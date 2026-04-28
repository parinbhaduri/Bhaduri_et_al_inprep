### Calculate damage metrics from model run data ###

# load environment and packages
import Pkg
Pkg.activate(".")
Pkg.instantiate()

#using Distributed
#addprocs(10, exeflags="--project=$(Base.active_project())")
using Base.Threads

using CSV, DataFrames
using Statistics 
using HDF5
using ProgressMeter
using StatsBase



include(joinpath(pwd(),"workflow","src","data_include.jl"))

out_dir = joinpath(@__DIR__,"data","shap_runs")
flood_year = 2013
haz_size = "Low"

#filename, pop_filename = filter(file -> startswith(file, string(flood_year)), readdir(out_dir))
pop_filename = "study_pop_share_data.h5"
#=
h5file = h5open(joinpath(out_dir,filename), "r")
#pop_dat = h5file["pop_data"]
price_dat = h5file["price_data"]
println("Price Dataset size: ", size(price_dat))
println(h5file["price_vars"][:])
=#
pop_file = h5open(joinpath(out_dir,pop_filename), "r")
pop_share_dat = pop_file["pop_share_data"]
println("Pop Share Dataset size: ", size(pop_share_dat))
println(pop_file["column_names"][:])

#Load property damage and exposure estimates
phil_damages = DataFrame(CSV.File(joinpath(dirname(pwd()), "philadelphia-data","flood_hazard", "data","phil_flood_dmg_ens.csv")))
phil_damages_no_unc = DataFrame(CSV.File(joinpath(dirname(pwd()), "philadelphia-data","flood_hazard", "data","phil_flood_dmg_no_unc.csv")))
phil_prop_exp = DataFrame(CSV.File(joinpath(dirname(pwd()), "philadelphia-data","flood_hazard", "data","phil_flood_prop_exp.csv")))

dmg_bgs = intersect(unique(phil_damages[phil_damages[!,"naccs_loss_$(flood_year)"] .> 0.0, :bg_id]), pop_file["GEOID"][:]) #Remove BGs not contained in in study area

## Read in damage estimates
function dat_to_dict(df,bg_list,col_name)
    dat_dict = Dict{Int, Vector{Float64}}()
    for bg in bg_list
        dat_array = zeros(3)
        for occ_cat in 1:3
            try
                dat_array[occ_cat] = first(df[
                    (df.bg_id .== bg) .&& 
                    (df.income_cat .== occ_cat), 
                    Symbol(col_name)
                ])
            catch
                continue
            end
        end
        dat_dict[bg] = dat_array
    end
    return dat_dict
end
println("Storing Exposed Property Counts (no uncertainty)...")
prop_exp_dict = dat_to_dict(phil_prop_exp, dmg_bgs,"num_prop_$(flood_year)")
println("Storing event damages (no uncertainty)...")
dam_no_unc_dict = dat_to_dict(phil_damages_no_unc, dmg_bgs,"naccs_loss_$(flood_year)")

println("Storing event damages (ensemble)...")
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
 

function process_damages(chunk, prop_exp_price, event_damages, event_dam_no_unc)
    batch_len = size(chunk, 1)
    
    # Vectorized house_pop (total people in each housing category)
    house_pop = zeros(batch_len, 3)
    for occ_cat in 1:3
        mask = chunk[:, :, 3] .== Float64(occ_cat)
        house_pop[:, occ_cat] = sum(chunk[:, :, 7] .* mask, dims=2)
    end
    
    # bg_pop calculation (total people in each income group (cols) in each housing category (rows) across all realizations)
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
    flood_damages = zeros(500, 3, batch_len)
    flood_dam_no_unc = zeros(3, batch_len)
    prop_exp_prices = zeros(3, batch_len)

    for j in 1:batch_len
        pop_prop = bg_pop[j, :, :] ./ house_pop[j, :]
        pop_prop[isnan.(pop_prop)] .= 0 #Remove nan if house_pop element is 0
        flood_damages[:, :, j] = event_damages * pop_prop
        flood_dam_no_unc[:, j] = event_dam_no_unc' * pop_prop
        prop_exp_prices[:, j] = prop_exp_price[j,:]' * pop_prop
    end
    
    return flood_damages, flood_dam_no_unc, prop_exp_prices
end

#findall(==(nothing),indexin(dmg_bgs,pop_file["GEOID"][:]))
### MAIN SCRIPT ###
n_areas = length(dmg_bgs)
batch_size = 1000
n_realizations = size(pop_share_dat,1)
n_batches = ceil(Int, n_realizations / batch_size)

# Create output HDF5 file
dat_dir = joinpath(out_dir,"post_process","flood_loss",haz_size)

mkpath(dat_dir)
loss_file = joinpath(dat_dir,"$(flood_year)_flood_loss.h5")

println("Starting Flood Loss Calculation for Flood Event $(flood_year)...")

h5open(loss_file, "w") do output_file
    
    # Pre-create dataset with chunking
    flood_loss = create_dataset(output_file, "flood_loss", 
                                     datatype(Float64), 
                                     dataspace(500, 3, n_realizations, n_areas),
                                     chunk=(500, 3, 1000, 1), deflate=9, shuffle=true
    )
    
    no_unc = create_dataset(output_file, "no_unc estimates", datatype(Float64), dataspace(n_realizations, n_areas, 6),
                        chunk=(1000,1,6), deflate=9, shuffle=true
    )
    # Metadata
    write(output_file, "loss columns", ["low income loss", "middle income loss","high income loss"])
    write(output_file, "no_unc columns", ["exposed value (low)", "exposed value (middle)","exposed value (high)","low income loss", "middle income loss","high income loss"])
    write(output_file, "loss_dims", ["damage realizations", "loss groups","ABM realizations", "block groups"])
    write(output_file, "no_unc_dims", ["ABM realizations", "block groups", "data columns"])
    write(output_file, "GEOID", dmg_bgs)
    
    #=
    # Write total housing prices
    ## Collect housing prices and household incomes within each damaged block group
    println("Starting housing price collection...")
    h5open(joinpath(out_dir,filename), "r") do price_file
        chunk_size = 6050
        total_rows = 159500
        price_dat = price_file["price_data"]
        for var_col in [1,2,3]
            for i in 1:chunk_size:total_rows
                end_idx = min(i + chunk_size - 1, total_rows)
                chunk = price_dat[i:end_idx, :, 5, var_col]
                cap_chunk = price_dat[i:end_idx, :, 5, var_col+3]
                for (j,ind) in enumerate(flpn_index)
                    bg_data[i:end_idx,j,var_col] = chunk[:,ind] .* cap_chunk[:,ind]
                end
                # Clear memory
                chunk = nothing
                cap_chunk = nothing
                GC.gc()
            end 
        end
    end
    println("Finished housing price collection!")
    =#
    # Open input file
    h5open(joinpath(out_dir,pop_filename), "r") do input_file
        pop_data = input_file["pop_share_data"]
        h5open(joinpath(out_dir,"post_process","flood_loss","study_bg_data.h5"), "r") do price_file
            price_data = price_file["bg_data"]
            # Loop over batches (read data once per batch)
            for batch_idx in 1:n_batches
                start_idx = (batch_idx - 1) * batch_size + 1
                end_idx = min(batch_idx * batch_size, n_realizations)
                batch_range = start_idx:end_idx
                
                println("Processing batch $(batch_idx)/$(n_batches) (realizations $(start_idx):$(end_idx))...")
                
                # Read this batch ONCE from HDF5
                chunk = pop_data[batch_range, :, :]
                
                price_chunk = price_data[batch_range, :, 1:3]
                
                # Process each area with this batch of data
                Threads.@threads for (i, bg) in collect(enumerate(dmg_bgs))
                    event_damages = event_damages_dict[bg]
                    event_dam_no_unc = dam_no_unc_dict[bg]
                    prop_exp_bg = prop_exp_dict[bg] 

                    matching_indices = findall(x->x == Float64(bg), chunk[:,:,1])
                    indices_by_row = [Int[] for _ in 1:size(chunk, 1)]
                    for idx in matching_indices
                        push!(indices_by_row[idx[1]], idx[2])
                    end

                    # Preallocate result
                    t_c = Array{eltype(chunk)}(undef, size(chunk, 1), 9, size(chunk,3))
                    
                    for b in eachindex(1:size(chunk,1))
                        t_c[b, :, :] = chunk[b, indices_by_row[b], :]
                    end

                    price_bg_ind = findfirst(==(bg), price_file["GEOID"][:])
                    t_p = price_chunk[:,price_bg_ind,1:3]

                    #Calculate total price value exposed 
                    prop_exp_price = t_p .* prop_exp_bg'

                    # Process this area for this batch
                    flood_loss_batch, flood_loss_batch_no_unc, prop_exp_prices = process_damages(t_c, prop_exp_price, event_damages, event_dam_no_unc)

                    flood_loss[:, :, batch_range, i] = flood_loss_batch
                    no_unc[batch_range, i, 1:3] = prop_exp_prices'
                    no_unc[batch_range, i, 4:6] = flood_loss_batch_no_unc'
                end
            end
            
            # Clear memory
            chunk = nothing
            GC.gc()
            
            #println("  Completed batch $(batch_idx)/$(n_batches)")
        end
    end
end

println("Done! Results saved to $(flood_year)_flood_loss.h5")

#=
chunk = reshape(pop_share_dat[49390, :, :], 1,size(pop_share_dat[49390, :, :])...)
tot_dam = zeros(500,3,1,length(dmg_bgs))
for (i,bg) in enumerate(dmg_bgs)
    matching_indices = findall(x->x == Float64(bg), chunk[:,:,1])
    indices_by_row = [Int[] for _ in 1:size(chunk, 1)]
    for idx in matching_indices
        push!(indices_by_row[idx[1]], idx[2])
    end

    # Preallocate result
    t_c = Array{eltype(chunk)}(undef, size(chunk, 1), 9, size(chunk,3))
                    
    for b in eachindex(1:size(chunk,1))
        t_c[b, :, :] = chunk[b, indices_by_row[b], :]
    end

    batch_len = size(chunk, 1)
        
    # Vectorized house_pop (total people in each housing category)
    house_pop = zeros(batch_len, 3)
    for occ_cat in 1:3
        mask = t_c[:, :, 3] .== Float64(occ_cat)
        house_pop[:, occ_cat] = sum(t_c[:, :, 7] .* mask, dims=2)
    end
        
    # bg_pop calculation (total people in each income group in each housing category across all realizations)
    bg_pop = zeros(batch_len, 3, 3)
    for occ_cat in 1:3
        for inc_group in 1:3
            mask = (t_c[:, :, 2] .== inc_group) .& (t_c[:, :, 3] .== occ_cat)
            for i in 1:batch_len
                idx = findfirst(mask[i, :])
                if !isnothing(idx)
                    bg_pop[i, occ_cat, inc_group] = t_c[i, idx, 7]
                end
            end
        end
    end
        
    # Calculate flood_loss for this batch
    flood_damages = zeros(500, 3, batch_len)
    for j in 1:batch_len
        pop_prop = bg_pop[j, :, :] ./ house_pop[j, :]
        if sum(pop_prop .> 1.0) >= 1.0
            println("block group id:", bg)
            println("block group order:", i)
            println(pop_prop)
        end
        flood_damages[:, :, j] = event_damages * pop_prop
    end

    tot_dam[:,:,:,i] = flood_damages
end

sum(filter(!isnan, tot_dam[445,:,1,:]))
tot_dam[445,:,1,:]
sum(phil_damages[phil_damages.sow_ind .== 445, :naccs_loss_2011])
=#
#=
bg_dat = h5open(joinpath(out_dir,"post_process","flood_loss","study_bg_data.h5"), "r")
price_data = bg_dat["bg_data"]

chunk = pop_share_dat[1:5, :, :]#reshape(pop_share_dat[49390, :, :], 1,size(pop_share_dat[49390, :, :])...)
price_chunk = price_data[1:5,:,1:3]#reshape(price_data[49390,:,1:3], 1,size(price_data[49390,:,1:3])...)
bg = dmg_bgs[1]

price_bg_ind = findfirst(==(bg), bg_dat["GEOID"][:])
t_p = price_chunk[:,price_bg_ind,1:3]

event_damages = event_damages_dict[bg]
prop_exp_bg = prop_exp_dict[bg]
event_dam_no_unc = dam_no_unc_dict[bg]
#Calculate total price value exposed 
prop_exp_price = t_p .* prop_exp_bg'

matching_indices = findall(x->x == Float64(bg), chunk[:,:,1])
indices_by_row = [Int[] for _ in 1:size(chunk, 1)]
for idx in matching_indices
    push!(indices_by_row[idx[1]], idx[2])
end

# Preallocate result
t_c = Array{eltype(chunk)}(undef, size(chunk, 1), 9, size(chunk,3))
                    
for b in eachindex(1:size(chunk,1))
    t_c[b, :, :] = chunk[b, indices_by_row[b], :]
end

batch_len = size(chunk, 1)


# Vectorized house_pop (total people in each housing category)
house_pop = zeros(batch_len, 3)
for occ_cat in 1:3
    mask = t_c[:, :, 3] .== Float64(occ_cat)
    house_pop[:, occ_cat] = sum(t_c[:, :, 7] .* mask, dims=2)
end
        
# bg_pop calculation (total people in each income group in each housing category across all realizations)
bg_pop = zeros(batch_len, 3, 3)
for occ_cat in 1:3
    for inc_group in 1:3
        mask = (t_c[:, :, 2] .== inc_group) .& (t_c[:, :, 3] .== occ_cat)
        for i in 1:batch_len
            idx = findfirst(mask[i, :])
            if !isnothing(idx)
                bg_pop[i, occ_cat, inc_group] = t_c[i, idx, 7]
            end
        end
    end
end
        
# Calculate flood_loss for this batch
flood_damages = zeros(500, 3, batch_len)
flood_dam_no_unc = zeros(3, batch_len)
prop_exp_prices = zeros(3, batch_len)

for j in 1:batch_len
    pop_prop = bg_pop[j, :, :] ./ house_pop[j, :]
    pop_prop[isnan.(pop_prop)] .= 0 #Remove nan if house_pop element is 0
    if sum(pop_prop .> 1.0) >= 1.0
        println("block group id:", bg)
        println("block group order:", i)
        println(pop_prop)
    end
    flood_damages[:, :, j] = event_damages * pop_prop
    flood_dam_no_unc[:, j] = event_dam_no_unc' * pop_prop
    #prop_exp_prices[:, j] = prop_exp_price[j,:]' * pop_prop
end
flood_damages[4,:,:]
flood_dam_no_unc

fl_batch, fl_no_unc, pr_prices = process_damages(t_c, prop_exp_price, event_damages, event_dam_no_unc)

fl_batch[4,:,:]
fl_no_unc
dam_no_unc[1:5,10,4:6]

no_unc =zeros(159500, 232, 6)
no_unc[1:5, 10, 4:6]# = fl_no_unc'
dam_dat[4,:,1:5,10]

pr_prices
dam_no_unc[1:5, 10, 1:3]
=#
### Calculate total flood losses in exposed areas

using Parquet2

function process_loss(dataset, price_dset; chunk_size=1000,total_rows=159500, var_col=1) #, index=flpn_index
    output = zeros(total_rows, size(dataset,1),2)
    for i in 1:chunk_size:total_rows
        end_idx = min(i + chunk_size - 1, total_rows)
        chunk = dataset[:, var_col, i:end_idx, :]
        price_chunk = price_dset[i:end_idx,:,var_col]
        
        #sum damages and prices across area
        flood_loss = dropdims(sum(x -> isnan(x) ? 0 : x, chunk, dims=3), dims=3)
        tot_price = dropdims(sum(x -> isnan(x) ? 0 : x, price_chunk, dims=2), dims=2)
        outcome = flood_loss ./ tot_price'
        replace!(outcome, NaN => 0) #NaN values occur when no agents exist in area. Assume burden is 0
        # Add chunk data to output
        output[i:end_idx,:,1] = flood_loss'
        output[i:end_idx,:,2] = outcome'
        # Clear memory
        chunk = nothing
        price_chunk = nothing
        GC.gc()
    end 
    return output
end

function convert_to_df(A)
    # Create index vectors
    dim1_indices = repeat(1:size(A, 1), outer=size(A, 2))
    dim2_indices = repeat(1:size(A, 2), inner=size(A, 1))

    df = DataFrame(
        model = dim1_indices,
        DDF = dim2_indices,
        loss_value = vec(A[:,:,1]),
        burden_value = vec(A[:,:,2])
    )
    return df
end

function spatial_tot_loss(array,price_dset,geo_list)
    spatial_summ = zeros(size(array,4),4)
    spatial_no_unc = zeros(size(array,4),2)
    for bg in axes(array,4)
        # Start w/uncertainty
        bg_loss = array[:,1:3,:,bg]
        tot_loss = dropdims(sum(x -> isnan(x) ? 0 : x, bg_loss, dims=2), dims=2)
        tot_mean = mean(skipmissing(tot_loss))
        tot_med = median(skipmissing(tot_loss))
        tot_quantiles = quantile(skipmissing(tot_loss), [0.025, 0.975])

        spatial_summ[bg,:] = [tot_mean,tot_med,tot_quantiles...]
        # w/o uncertainty
        bg_loss_no_unc = price_dset[:,bg,4:6]
        tot_loss_no_unc = dropdims(sum(x -> isnan(x) ? 0 : x, bg_loss_no_unc, dims=2), dims=2)
        tot_mean_no_unc = mean(skipmissing(tot_loss_no_unc))
        tot_med_no_unc = median(skipmissing(tot_loss_no_unc))

        spatial_no_unc[bg,:] = [tot_mean_no_unc, tot_med_no_unc]
        #=
        bg_burden = bg_loss ./ bg_price'
        replace!(bg_burden, NaN => 0) #NaN values occur when no agents exist in bg. Assume burden is 0
        burd_med = median(skipmissing(bg_burden))
        burd_quantiles = quantile(skipmissing(bg_burden), [0.025, 0.975])
        spatial_summ[bg,:,2] = [burd_med,burd_quantiles...]
        =#
    end
    spatial_loss = DataFrame(GEOID = geo_list, loss_mean = spatial_summ[:,1,1], loss_med = spatial_summ[:,2],
                        loss_lower = spatial_summ[:,3], loss_upper = spatial_summ[:,4],
                        loss_mean_no_unc =  spatial_no_unc[:,1], loss_med_no_unc =  spatial_no_unc[:, 2]
    )
    return spatial_loss
end

Haz_Dict = Dict(
    "High" => [2011],#,1989,1996], #,
    "Medium" => [2018]#1991,2018,1981],
    #"Low" => [1988,2013,2010],
)

for haz_size in ["Medium"]#"High",=, "Low"]
    for flood_year in Haz_Dict[haz_size]
        dat_dir = joinpath(out_dir,"post_process","flood_loss",haz_size)
        dam_file = h5open(joinpath(dat_dir,"$(flood_year)_flood_loss.h5"), "r")
        dam_dat = dam_file["flood_loss"]
        dam_no_unc = dam_file["no_unc estimates"]
        #=
        #bg_dat = dam_file["bg_data"]
        #println("Damage Dataset size: ", size(dam_dat))
        #println("No Uncertainty Estimate Dataset size: ", size(dam_no_unc))

        println("Starting Loss Estimates for $(flood_year)...")
        flpn_low_loss = process_loss(dam_dat, dam_no_unc; var_col=1)
        #flpn_low_df = convert_to_df(flpn_low_loss)
        #Parquet2.writefile(joinpath(dat_dir,"$(flood_year)_model_outcome_flpn_loss_low.parquet"), flpn_low_df)
        #CSV.write(joinpath(dat_dir,"$(flood_year)_model_outcome_flpn_loss_low.csv"), flpn_low_df)
        flpn_med_loss = process_loss(dam_dat, dam_no_unc; var_col=2)
        #flpn_med_df = convert_to_df(flpn_med_loss)
        #Parquet2.writefile(joinpath(dat_dir,"$(flood_year)_model_outcome_flpn_loss_med.parquet"), flpn_med_df)
        flpn_high_loss = process_loss(dam_dat, dam_no_unc; var_col=3)
        #flpn_high_df = convert_to_df(flpn_high_loss)
        #Parquet2.writefile(joinpath(dat_dir,"$(flood_year)_model_outcome_flpn_loss_high.parquet"), flpn_high_df)

        flpn_loss_tot = flpn_low_loss .+ flpn_med_loss .+ flpn_high_loss
        #Fix burden estimate
        for i in 1:1000:159500
            end_idx = min(i + 1000 - 1, 159500)
            price_chunk = dam_no_unc[i:end_idx,:,1:3] #
            tot_price = dropdims(sum(x -> isnan(x) ? 0 : x, price_chunk, dims=(2,3)), dims=(2,3))

            outcome = flpn_loss_tot[i:end_idx,:,1] ./ tot_price #i:end_idx
            replace!(outcome, NaN => 0) #NaN values occur when no agents exist in area. Assume burden is 0
            flpn_loss_tot[i:end_idx,:,2] = outcome
        end
        flpn_df = convert_to_df(flpn_loss_tot)
        Parquet2.writefile(joinpath(dat_dir,"$(flood_year)_model_outcome_flpn_loss.parquet"), flpn_df)
        =#
        #Calculate total losses by block group 
        println("Starting Spatial Loss Estimates for $(flood_year)...")
        total_spatial_loss = spatial_tot_loss(dam_dat, dam_no_unc,dam_file["GEOID"][:])
        CSV.write(joinpath(dat_dir,"$(flood_year)_model_outcome_spatial_loss.csv"),total_spatial_loss)
    end
end

#=

## Calculate Concentration curves of flood losses by income
study_bg_file = h5open(joinpath(out_dir,"post_process","flood_loss","study_bg_data.h5"), "r")
study_bg_dat = study_bg_file["bg_data"]

sort_inc_file = h5open(joinpath(out_dir,"post_process","flood_loss","sorted_inc_indices.h5"),"r")
all_study_ind = sort_inc_file["bg_indices"][:,:]
dmg_study_ind = sort_inc_file["damaged_indices"][:,:]
inc_ind = sort_inc_file["inc_indices"][:,:]

# Create reverse mapping: All bgs to dmg bgs for event only
sub_indices = filter(!isnothing,indexin(dam_file["GEOID"][:], sort_inc_file["Study Area BGs"][:]))
bg_map = Dict{Int, Union{Int, Nothing}}()
for (new_idx, old_idx) in enumerate(sub_indices)
    bg_map[old_idx] = new_idx
end
sort_inc_file["Study Area BGs"][36]
sort_inc_file["Damaged BGs"][3]
# Create reverse mapping: All damage bgs to dmg bgs for event only
sub_dam_indices = filter(!isnothing,indexin(dam_file["GEOID"][:], sort_inc_file["Damaged BGs"][:]))
dam_bg_map = Dict{Int, Union{Int, Nothing}}()
for (new_idx, old_idx) in enumerate(sub_dam_indices)
    dam_bg_map[old_idx] = new_idx
end
l = 5
m_r = 1
f_a = all_study_ind[m_r,l]
f_d = dmg_study_ind[m_r,l]
f_i = inc_ind[m_r,l]
# Collect damage values
damages = dam_dat[:,f_i,m_r,dam_bg_map[f_d]]
# Collect housing price value
price = dam_no_unc[m_r,dam_bg_map[f_d], f_i]
# Collect population values
pop = study_bg_dat[m_r, bg_map[f_a], f_i + 6]

#Calculate Risk burden
burden = damages ./ price
#Calculate Summary Statistics

for i in 1:159500
    # Filter and remap in one pass
    j_200_row = Int[]
    k_200_row = Int[]
        
    for col in 1:696
        j_val = j_232[i, col]
        if haskey(index_map_232_to_200, j_val)
            push!(j_200_row, index_map_232_to_200[j_val])
            push!(k_200_row, k[i, col])
        end
    end
end


#access loss data
loss_sl = zeros(500,936)
dam_sl = dam_dat[:, :, 1, :]
for i in eachindex(location_indices[1, :])
    loss_sl[:, i] = dam_sl[:, data_indices[1, i], location_indices[1, i]]
end



for i in 1:159500
    # Extract only the selected indices from the i-th slice
    slice = data[i, selected_indices, :]  # Now 312 x 6
        
    # Extract columns 4-6 for sorting
    sort_keys = slice[:, 4:6]
        
    # Get the sorting permutation based on columns 4, 5, 6
    perm = sortperm([sort_keys[j, :] for j in 1:312])
        
    # Store the permutation (these are indices into selected_indices)
    sort_orders[i] = perm
        
    # Sort the slice according to this permutation
    sorted_slice = slice[perm, :]
        
    # Flatten to 1D (312 * 6 = 1872 elements)
    result[i, :] = vec(sorted_slice')
end
=#

## Calculate Median flood loss and 95% uncertainty interval for each block group
function spatial_loss(array,price_dset,geo_list;var_col=1)
    spatial_summ = zeros(size(array,4),3,2)
    for bg in axes(array,4)
        bg_loss = array[:, var_col, :, bg]
        bg_price = price_dset[:,bg,var_col]
        bg_med = median(skipmissing(bg_loss))
        bg_quantiles = quantile(skipmissing(bg_loss), [0.025, 0.975])
        spatial_summ[bg,:,1] = [bg_med,bg_quantiles...]

        bg_burden = bg_loss ./ bg_price'
        replace!(bg_burden, NaN => 0) #NaN values occur when no agents exist in bg. Assume burden is 0
        burd_med = median(skipmissing(bg_burden))
        burd_quantiles = quantile(skipmissing(bg_burden), [0.025, 0.975])
        spatial_summ[bg,:,2] = [burd_med,burd_quantiles...]
    end
    spatial_loss = DataFrame(GEOID = geo_list, loss_med = spatial_summ[:,1,1], burden_med = spatial_summ[:,1,2],
                        loss_lower = spatial_summ[:,2,1], burden_lower = spatial_summ[:,2,2],
                        loss_upper = spatial_summ[:,3,1], burden_upper = spatial_summ[:,3,2]
    )
    return spatial_loss
end

println("Calculating spatial uncertainty in flood loss for low income...")
sp_loss_low = spatial_loss(dam_dat,dam_no_unc,dam_file["GEOID"][:];var_col=1)
CSV.write(joinpath(dat_dir,"$(flood_year)_model_outcome_spatial_loss_low.csv"),sp_loss_low)
println("Calculating spatial uncertainty in flood loss for middle income...")
sp_loss_med = spatial_loss(dam_dat,dam_no_unc,dam_file["GEOID"][:];var_col=2)
CSV.write(joinpath(dat_dir,"$(flood_year)_model_outcome_spatial_loss_med.csv"),sp_loss_med)
println("Calculating spatial uncertainty in flood loss for high income...")
sp_loss_high = spatial_loss(dam_dat,dam_no_unc,dam_file["GEOID"][:];var_col=3)
CSV.write(joinpath(dat_dir,"$(flood_year)_model_outcome_spatial_loss_high.csv"),sp_loss_high)


dat_dir = joinpath(out_dir,"post_process","flood_loss","High")
dam_file = h5open(joinpath(dat_dir,"2011_flood_loss.h5"), "r")
dam_dat = dam_file["flood_loss"]
dam_no_unc = dam_file["no_unc estimates"]
## Repeat for Total flood losses w/ and w/o uncertainty
function spatial_tot_loss(array,price_dset,geo_list)
    spatial_summ = zeros(size(array,4),4)
    spatial_no_unc = zeros(size(array,4),2)
    for bg in axes(array,4)
        # Start w/uncertainty
        bg_loss = array[:,1:3,:,bg]
        tot_loss = dropdims(sum(x -> isnan(x) ? 0 : x, bg_loss, dims=2), dims=2)
        tot_mean = mean(skipmissing(tot_loss))
        tot_med = median(skipmissing(tot_loss))
        tot_quantiles = quantile(skipmissing(tot_loss), [0.025, 0.975])

        spatial_summ[bg,:] = [tot_mean,tot_med,tot_quantiles...]
        # w/o uncertainty
        bg_loss_no_unc = price_dset[:,1,4:6]
        tot_loss_no_unc = dropdims(sum(x -> isnan(x) ? 0 : x, bg_loss_no_unc, dims=2), dims=2)
        tot_mean_no_unc = mean(skipmissing(tot_loss_no_unc))
        tot_med_no_unc = median(skipmissing(tot_loss_no_unc))

        spatial_no_unc[bg,:] = [tot_mean_no_unc, tot_med_no_unc]
        #=
        bg_burden = bg_loss ./ bg_price'
        replace!(bg_burden, NaN => 0) #NaN values occur when no agents exist in bg. Assume burden is 0
        burd_med = median(skipmissing(bg_burden))
        burd_quantiles = quantile(skipmissing(bg_burden), [0.025, 0.975])
        spatial_summ[bg,:,2] = [burd_med,burd_quantiles...]
        =#
    end
    spatial_loss = DataFrame(GEOID = geo_list, loss_mean = spatial_summ[:,1,1], loss_med = spatial_summ[:,2],
                        loss_lower = spatial_summ[:,3], loss_upper = spatial_summ[:,4],
                        loss_mean_no_unc =  spatial_no_unc[:,1], loss_mean_no_unc =  spatial_no_unc[:, 2]
    )
    return spatial_loss
end
bg_loss = dam_dat[:,1:3,:,1]
tot_loss = dropdims(sum(x -> isnan(x) ? 0 : x, bg_loss, dims=2), dims=2)
tot_mean = mean(skipmissing(tot_loss))
tot_med = median(skipmissing(tot_loss))
tot_quantiles = quantile(skipmissing(tot_loss), [0.025, 0.975])

#No uncertainty
bg_loss_no_unc = dam_no_unc[:,1,4:6]
tot_loss_no_unc = dropdims(sum(x -> isnan(x) ? 0 : x, bg_loss_no_unc, dims=2), dims=2)
tot_mean_no_unc = mean(skipmissing(tot_loss_no_unc))
tot_med_no_unc = median(skipmissing(tot_loss_no_unc))
#tot_quantiles_no_unc = quantile(skipmissing(tot_loss_no_unc), [0.025, 0.975])
#=
b_l = dropdims(sum(x -> isnan(x) ? 0 : x, dam_dat[:, 1:3, :, 1], dims=2), dims=2)
b_u = dropdims(sum(x -> isnan(x) ? 0 : x, dam_no_unc[:,1,4:6], dims=2), dims=2)
#Calculate % diff
dam_diff = (b_l .- b_u') ./ b_u'
diff_med = median(skipmissing(dam_diff))
diff_quantiles = quantile(skipmissing(dam_diff), [0.025, 0.975])
function loss_diff(array,no_unc_array,geo_list;var_col=1)
    spatial_summ = zeros(size(array,4),3)
    for bg in axes(array,4)
        if bg % 10 == 0
            println("starting bg no. $(bg)...")
        end
        bg_loss = dropdims(sum(x -> isnan(x) ? 0 : x, array[:, 1:3, :, bg], dims=2), dims=2) #total loss
        bg_loss_no_unc = dropdims(sum(x -> isnan(x) ? 0 : x, no_unc_array[:,bg,4:6], dims=2), dims=2) # total loss no_unc
        #Calculate % diff
        dam_diff = (bg_loss .- bg_loss_no_unc') ./ bg_loss_no_unc'
        replace!(dam_diff, NaN => 0)
        diff_med = median(skipmissing(dam_diff))
        diff_quantiles = quantile(skipmissing(dam_diff), [0.025, 0.975])
        spatial_summ[bg,:] = [diff_med,diff_quantiles...]
    end
    spatial_loss = DataFrame(GEOID = geo_list, dam_diff_med = spatial_summ[:,1,1], 
                        dam_diff_lower = spatial_summ[:,2],
                        dam_diff_upper = spatial_summ[:,3], 
    )
    return spatial_loss
end
diff_loss = loss_diff(dam_dat,dam_no_unc,dam_file["GEOID"][:];var_col=1)
CSV.write(joinpath(dat_dir,"$(flood_year)_model_outcome_spatial_loss_diff.csv"),diff_loss)
=#
println("data collection finished!")

