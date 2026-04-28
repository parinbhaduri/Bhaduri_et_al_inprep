### File creates a hdf5 file storing housing price and household income data
# across all realizations during the time of the flood shock ###
# load environment and packages
import Pkg
Pkg.activate(".")
Pkg.instantiate()

using CSV, DataFrames
using Statistics 
using HDF5
using ProgressMeter
using StatsBase

using Base.Threads

include(joinpath(pwd(),"workflow","src","data_include.jl"))

out_dir = joinpath(@__DIR__,"data","shap_runs")
flood_year = 2011
filename = first(filter(file -> startswith(file, (string(flood_year)*"_abm")), readdir(out_dir)))
pop_filename = "study_pop_share_data.h5"

h5file = h5open(joinpath(out_dir,filename), "r")
#pop_dat = h5file["pop_data"]
price_dat = h5file["price_data"]
println("Price Dataset size: ", size(price_dat))
println(h5file["price_vars"][:])

pop_file = h5open(joinpath(out_dir,pop_filename), "r")
pop_share_dat = pop_file["pop_share_data"]
println("Pop Share Dataset size: ", size(pop_share_dat))
println(pop_file["column_names"][:])


## MAIN SCRIPT
study_bgs = pop_file["GEOID"][:]
#matching_indices_dict = Dict(bg => findall(x->x == Float64(bg), GEOIDs) for bg in study_bgs)

n_areas = length(study_bgs)
batch_size = 1000
n_realizations = size(pop_share_dat,1)
n_batches = ceil(Int, n_realizations / batch_size)
# Create output HDF5 file
dat_dir = joinpath(out_dir,"post_process","flood_loss")
bg_file = joinpath(dat_dir,"study_bg_data.h5")

h5open(bg_file, "w") do output_file
    
    # Pre-create dataset with chunking
    bg_data = create_dataset(output_file, "bg_data", datatype(Float64), dataspace(n_realizations, n_areas, 9),
                        chunk=(1000,1,9), deflate=9, shuffle=true
    )
    # Metadata
    write(output_file, "block group columns", ["low house price", "middle house price","high house price", "low income", "middle income","high income","low inc. pop", "middle inc. pop","high inc. pop"])
    write(output_file, "dimensions", ["model realizations", "block groups","data columns"])
    write(output_file, "GEOID", study_bgs)
    
    # Write total housing prices
    println("Starting housing price collection...")
    h5open(joinpath(out_dir,filename), "r") do price_file
        chunk_size = 6050
        total_rows = 159500
        price_dat = price_file["price_data"]


        for var_col in [1,2,3]
            println("  Processing price variable $var_col/3...")
            for i in 1:chunk_size:total_rows
                end_idx = min(i + chunk_size - 1, total_rows)
                
                # Read and compute in one go (Only need avg_price, but uncomment lines if total housing value is needed)
                bg_data[i:end_idx, :, var_col] = price_dat[i:end_idx, :, 5, var_col] #.* 
                                                #price_dat[i:end_idx, :, 5, var_col+3]
            end 
            #for i in 1:chunk_size:total_rows
                #end_idx = min(i + chunk_size - 1, total_rows)
    
                #chunk = price_dat[i:end_idx, :, 5, var_col] #First column is initial state, so year 4 is column 5
                #cap_chunk = price_dat[i:end_idx, :, 5, var_col+3]
                
                #for ind in 1:length(study_bgs)
                #    bg_data[i:end_idx,ind,var_col] = chunk[:,ind] .* cap_chunk[:,ind]
                #end
            #end 
        end
    end
    println("Finished housing price collection!")
    
    # Write total household income
    println("Starting household income collection with $(Threads.nthreads()) threads...")
    h5open(joinpath(out_dir,pop_filename), "r") do pop_file
        pop_dat = pop_file["pop_share_data"]
        # Create one set of buffers PER THREAD to avoid contention
        thread_buffers = [
            (income_ag_pop = zeros(Float64, batch_size, 3), #total number of agents in income group
            income_pop = zeros(Float64, batch_size, 3), #total number of people in income group
            income = zeros(Float64, batch_size, 3), #total income in housing group
            avg_income = zeros(Float64, batch_size, 3))
            for _ in 1:Threads.nthreads()
        ]
        # Loop over batches (read data once per batch)
        for batch_idx in 1:n_batches
            start_idx = (batch_idx - 1) * batch_size + 1
            end_idx = min(batch_idx * batch_size, n_realizations)
            batch_range = start_idx:end_idx
            batch_len = length(batch_range)
            
            println("Processing batch $(batch_idx)/$(n_batches) (realizations $(start_idx):$(end_idx))...")
            # Read this batch ONCE from HDF5
            chunk = pop_dat[batch_range, :, :]

            Threads.@threads for i in eachindex(study_bgs)
                # Get this thread's ID and its buffers
                tid = Threads.threadid()
                buffers = thread_buffers[tid]
                
                bg = study_bgs[i]
                matching_indices = findall(x->x == Float64(bg), chunk[:,:,1])
                indices_by_row = [Int[] for _ in 1:size(chunk, 1)]
                for idx in matching_indices
                    push!(indices_by_row[idx[1]], idx[2])
                end

                # Preallocate result
                chunk_matches = Array{eltype(chunk)}(undef, size(chunk, 1), 9, size(chunk,3))
                
                for b in eachindex(1:size(chunk,1))
                    chunk_matches[b, :, :] = chunk[b, indices_by_row[b], :]
                end

                # Reset this thread's buffers
                fill!(view(buffers.income_ag_pop, 1:batch_len, :), 0.0)
                fill!(view(buffers.income_pop, 1:batch_len, :), 0.0)
                fill!(view(buffers.income, 1:batch_len, :), 0.0)
             
                for inc_group in 1:3
                    @views begin
                        mask = chunk_matches[:, :, 2] .== Float64(inc_group)
                        buffers.income_ag_pop[1:batch_len, inc_group] = sum(chunk_matches[:, :, 4] .* mask, dims=2)
                        buffers.income[1:batch_len, inc_group] = sum(chunk_matches[:, :, 5] .* mask, dims=2)

                        buffers.income_pop[1:batch_len, inc_group] = sum(chunk_matches[:, :, 7] .* mask, dims=2)
                    end
                end
                #Safe Division in case no people present in block group
                @views for inc_group in 1:3 
                    for j in 1:batch_len
                        if buffers.income_ag_pop[j, inc_group] > 0
                            buffers.avg_income[j, inc_group] = buffers.income[j, inc_group] / buffers.income_ag_pop[j, inc_group]
                        else
                            buffers.avg_income[j, inc_group] = 0.0
                        end
                    end
                end

                bg_data[batch_range, i, 4:6] = view(buffers.avg_income, 1:batch_len, :)
                bg_data[batch_range, i, 7:9] = view(buffers.income_pop, 1:batch_len, :)

            end

            # Clear memory
            chunk = nothing
            GC.gc()
        end
    end
    println("Finished household income collection!")
    
end

## Calculate Summary Demographic Characteristics for each block group


#=
###IGNORE ###
## Create index ordering of block groups from smallest to largest income by model realization 
println("Creating Sorted Indices by Income for Block Group Data Subsetting...")
#Create file to store ordering 
sort_inc_file = joinpath(dat_dir,"sorted_inc_indices.h5")
#Read in damages
phil_damages = DataFrame(CSV.File(joinpath(dirname(pwd()), "philadelphia-data","flood_hazard", "data","phil_flood_dmg_ens.csv")))
#Select block group indices that experience damages
dam_GEOID = unique(phil_damages.bg_id)


h5open(bg_file, "r") do infile
    bg_dat = infile["bg_data"]
    exp_indices = filter(!isnothing,indexin(dam_GEOID, infile["GEOID"][:]))
    exp_ind_size = length(exp_indices)
    flat_dim_size = exp_ind_size * 3

    h5open(sort_inc_file, "w") do outfile
        # Create datasets for indices
        bg_orig_dset = create_dataset(outfile, "bg_indices", datatype(Int32), dataspace(159500, flat_dim_size))  # For study area arrays (755 bgs)
        bg_subset_dset = create_dataset(outfile, "damaged_indices", datatype(Int32), dataspace(159500, flat_dim_size))  # For damage arrays (242 total bgs)
        inc_dset = create_dataset(outfile, "inc_indices", datatype(Int32), dataspace(159500, flat_dim_size)) #Income Group Columns 

        # Metadata
        write(outfile, "Study Area BGs", infile["GEOID"][:])
        write(outfile, "Damaged BGs", dam_GEOID)
        
        chunk_size = 1000
        
        for start_idx in 1:chunk_size:159500
            end_idx = min(start_idx + chunk_size - 1, 159500)
            chunk_rows = end_idx - start_idx + 1
            t_d = bg_dat[start_idx:end_idx,:,4:6]

            exp_chunk = t_d[:,exp_indices,:]
            reshaped_chunk = reshape(exp_chunk, chunk_rows, flat_dim_size)

            bg_orig_chunk = Matrix{Int}(undef, chunk_rows, flat_dim_size)
            bg_subset_chunk = Matrix{Int}(undef, chunk_rows, flat_dim_size)
            inc_chunk = Matrix{Int}(undef, chunk_rows, flat_dim_size)
            for i in 1:chunk_rows
                sorted_idx = sortperm(reshaped_chunk[i, :])
                            
                for col in 1:flat_dim_size
                    linear_idx = sorted_idx[col]
                    bg = mod1(linear_idx, exp_ind_size)
                    inc = div(linear_idx - 1, exp_ind_size) + 1
                    
                    # Store both: original 755-index and subset 232-index
                    bg_orig_chunk[i, col] = exp_indices[bg]
                    bg_subset_chunk[i, col] = bg
                    inc_chunk[i, col] = inc
                end
            end

            # Write chunks to file
            bg_orig_dset[start_idx:end_idx, :] = bg_orig_chunk
            bg_subset_dset[start_idx:end_idx, :] = bg_subset_chunk
            inc_dset[start_idx:end_idx, :] = inc_chunk
            
            println("Processed rows $(start_idx) to $(end_idx)")
        end
    end
end

println("Done!")
pop_dat = pop_file["pop_share_data"]
chunk = pop_dat[1:2, :, :]
batch_len = size(chunk, 1)
bg = study_bgs[300]
matching_indices = findall(x->x == Float64(bg), chunk[:,:,1])
indices_by_row = [Int[] for _ in 1:size(chunk, 1)]
for idx in matching_indices
    push!(indices_by_row[idx[1]], idx[2])
end
# Preallocate result
chunk_matches = Array{eltype(chunk)}(undef, size(chunk, 1), 9, size(chunk,3))
                
for b in eachindex(1:size(chunk,1))
    chunk_matches[b, :, :] = chunk[b, indices_by_row[b], :]
end

income_pop = zeros(Float64, 2, 3) #total number of agents in income group
income = zeros(Float64, 2, 3) #total income in housing group
avg_income = zeros(Float64, 2, 3)

for inc_group in 1:3
    @views begin
        mask = chunk_matches[:, :, 2] .== Float64(inc_group)
        income_pop[1:batch_len, inc_group] = sum(chunk_matches[:, :, 4] .* mask, dims=2)
        income[1:batch_len, inc_group] = sum(chunk_matches[:, :, 5] .* mask, dims=2)
    end
end
#Safe Division in case no people present in block group
@views for inc_group in 1:3 
    for j in 1:batch_len
        if income_pop[j, inc_group] > 0
            avg_income[j, inc_group] = income[j, inc_group] / income_pop[j, inc_group]
        else
            avg_income[j, inc_group] = 0.0
        end
    end
end
=#