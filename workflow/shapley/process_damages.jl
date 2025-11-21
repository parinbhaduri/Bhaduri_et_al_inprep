### Calculate damage metrics from model run data ###

# load environment and packages
import Pkg
Pkg.activate(".")
Pkg.instantiate()

using Distributed
addprocs(8, exeflags="--project=$(Base.active_project())")

@everywhere begin
    using CSV, DataFrames
    using Statistics 
    using HDF5
    using ProgressMeter
    using StatsBase
end


include(joinpath(pwd(),"workflow","src","data_include.jl"))

out_dir = joinpath(@__DIR__,"data","shap_runs")
filename = "2011_abm_data_142772.h5"
pop_filename = "1989_pop_share_data_DESKTOP.h5"

h5file = h5open(joinpath(out_dir,filename), "r")
#pop_dat = h5file["pop_data"]
price_dat = h5file["price_data"]
println("Price Dataset size: ", size(price_dat))
println(h5file["price_vars"][:])

pop_file = h5open(joinpath(out_dir,pop_filename), "r")
pop_share_dat = pop_file["pop_share_data"]
println("Pop Share Dataset size: ", size(pop_share_dat))
println(pop_file["column_names"][:])

#pop_share_data = pop_share_dat[:,:,:]
#close(pop_file)

#Subset Area to Floodplain (Select Block Groups)
phil_damages = DataFrame(CSV.File(joinpath(dirname(pwd()), "philadelphia-data","flood_hazard", "data","phil_flood_dmg_ens.csv")))
flood_year = read(h5file["historical flood year"])
exp_bgs = phil_flood_record[phil_flood_record[!,string(flood_year)] .> 0,"GEOID"]
dmg_bgs = intersect(exp_bgs,unique(phil_damages.bg_id))


# Read in damage estimates
println("Storing event damages...")
event_damages_dict = Dict{Int, Matrix{Float64}}()
for bg in dmg_bgs
    event_damages = zeros(500, 3)
    for occ_cat in 1:3
        try
            event_damages[:, occ_cat] = phil_damages[
                (phil_damages.bg_id .== bg) .&& 
                (phil_damages.income_cat .== occ_cat), 
                Symbol("naccs_loss_$(flood_year)")
            ]
        catch
            continue
        end
    end
    event_damages_dict[bg] = event_damages
end
#Calculate damage burden for each housing category in each block group
#(Damage burden is the proportion of total damages to total housing value)
   
#h5open(h5_filepath, "r") do file
function process_loss(chunk, event_damages)
    batch_len = size(chunk, 1)
    
    # Vectorized house_pop
    house_pop = zeros(batch_len, 3)
    for occ_cat in 1:3
        mask = chunk[:, :, 3] .== Float64(occ_cat)
        house_pop[:, occ_cat] = sum(chunk[:, :, 7] .* mask, dims=2)
    end
    
    # bg_pop calculation
    bg_pop = zeros(batch_len, 3, 3)
    for occ_cat in 1:3
        for inc_group in 1:3
            mask = (chunk[:, :, 2] .== inc_group) .& (chunk[:, :, 3] .== occ_cat)
            for i in 1:batch_len
                idx = findfirst(mask[i, :])
                if !isnothing(idx)
                    bg_pop[i, occ_cat, inc_group] = chunk[i, idx, 7]
                end
            end
        end
    end
    
    # Calculate flood_loss for this batch
    flood_loss = zeros(500, 3, batch_len)
    for i in 1:batch_len
        pop_prop = bg_pop[i, :, :] ./ house_pop[i, :]'
        flood_loss[:, :, i] = event_damages * pop_prop
    end
    
    return flood_loss
end

### MAIN SCRIPT ###
batch_size = 1000
n_realizations = size(pop_share_dat,1)
n_batches = ceil(Int, n_realizations / batch_size)
# Create output HDF5 file
h5open("flood_loss_results.h5", "w") do output_file
    # Pre-create dataset with chunking
    flood_loss_all = create_dataset(output_file, "flood_loss", 
                                     datatype(Float64), 
                                     dataspace(500, 3, n_realizations, n_areas),
                                     chunk=(500, 3, 1000, 1))
    
    # Open input file
    h5open(h5_filepath, "r") do input_file
        pop_data = input_file[dataset_name]
        
        # Loop over batches (read data once per batch)
        for batch_idx in 1:n_batches
            start_idx = (batch_idx - 1) * batch_size + 1
            end_idx = min(batch_idx * batch_size, n_realizations)
            batch_range = start_idx:end_idx
            
            println("Processing batch $(batch_idx)/$(n_batches) (realizations $(start_idx):$(end_idx))...")
            
            # Read this batch ONCE from HDF5
            chunk = pop_data[batch_range, :, :]
            
            # Process each area with this batch of data
            for area_idx in 1:232
                matching_indices = matching_indices_dict[area_idx]
                event_damages = event_damages_dict[area_idx]
                
                # Extract relevant columns for this area
                t_c = @view chunk[:, matching_indices, :]
                
                # Process this area for this batch
                flood_loss_batch = process_area_for_batch(t_c, event_damages)
                
                # Write results directly to disk
                flood_loss_all[:, :, batch_range, area_idx] = flood_loss_batch
            end
            
            println("  Completed batch $(batch_idx)/$(n_batches)")
        end
    end
end

println("Done! Results saved to flood_loss_results.h5")


















flood_loss = zeros(500, 3, 20, n_areas)
#for area_idx in n_areas
matching_indices = findall(x->x == Float64(dmg_bgs[1]), GEOIDs)

chunk = pop_share_dat[1:11, :, :]
t_c = chunk[:, matching_indices, :]
event_damages = zeros(500, 3)
for occ_cat in 1:3
    try
        event_damages[:, occ_cat] = phil_damages[
            (phil_damages.bg_id .== dmg_bgs[1]) .&& 
            (phil_damages.income_cat .== occ_cat), 
            Symbol("naccs_loss_1989")
        ]
    catch
        continue
    end
end
#calculate total pop for each housing category
for i in 1:11
    house_pop = zeros(3)
    for occ_cat in 1:3
        mask = t_c[i, :, 3] .== Float64(occ_cat)
        house_pop[occ_cat] = sum(t_c[i, :, 7] .* mask)
    end
    # Create 3x3 matrix of proportion of people in income category living in housing category
    #(rows=housing, col=income) 
    bg_pop = zeros(3, 3)
    for occ_cat in 1:3
        for inc_group in 1:3
            ind = findall(
                (t_c[i,:,2] .== inc_group) .&
                (t_c[i,:,3] .== occ_cat)
            )
            if !isempty(ind)
                bg_pop[occ_cat, inc_group] = t_c[i, ind[1], 7]
            end
        end
    end

    # Calculate pop_prop for this row
    pop_prop = bg_pop ./ house_pop'

    # Calculate and save flood_loss
    flood_loss[:, :, i, area_idx] = event_damages * pop_prop
end