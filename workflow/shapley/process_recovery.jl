### Calculate output metrics from model run data ###

# load environment and packages
import Pkg
Pkg.activate(".")
Pkg.instantiate()


using CSV, DataFrames
using Statistics 
using HDF5
using ProgressMeter
using StatsBase

include(joinpath(pwd(),"workflow","src","data_include.jl"))

function process_pop(dataset;chunk_size=6050,total_rows=159500, var_col = 1, subset=false, index=flpn_index)
    output = zeros(total_rows,size(dataset, 3))
    for i in 1:chunk_size:total_rows
        end_idx = min(i + chunk_size - 1, total_rows)
        chunk = dataset[i:end_idx, :, :, var_col]
        if subset
            outcome = zeros(size(chunk, 1),size(chunk, 3))
            for ind in index
                outcome += chunk[:,ind,:]
            end
        else
            outcome = dropdims(sum(chunk,dims=2), dims=2)
        end
        
        # Add chunk data to output
        output[i:end_idx, :] = outcome
        # Clear memory
        chunk = nothing
        GC.gc()
    end 
    return output
end

function process_price(dataset;chunk_size=6050,total_rows=159500, var_col = 1, subset=false, index=flpn_index)
    output = zeros(total_rows,size(dataset, 3))
    for i in 1:chunk_size:total_rows
        end_idx = min(i + chunk_size - 1, total_rows)
        chunk = dataset[i:end_idx, :, :, var_col]
        cap_chunk = dataset[i:end_idx, :, :, var_col+3]
        if subset
            outcome = zeros(size(chunk, 1),size(chunk, 3))
            total_cap = zeros(size(chunk, 1),size(chunk, 3))
            for ind in index
                outcome += chunk[:,ind,:] #.* cap_chunk[:,ind,:]
                total_cap += cap_chunk[:,ind,:]
            end
        else
            outcome = dropdims(sum(chunk,dims=2), dims=2) ./ dropdims(sum(cap_chunk,dims=2), dims=2)
        end
        
        # Add chunk data to output
        output[i:end_idx, :] = outcome ./ length(index)#total_cap
        # Clear memory
        chunk = nothing
        GC.gc()
    end 
    return output
end

function normal_set(dset)
    init_cond = dset[:,1]
    return ((dset .- init_cond) ./ init_cond) .* 100
end


###Read in data###
## Create dict to categorize flood years by hazard size
haz_cat = DataFrame(CSV.File(joinpath(dirname(pwd()), "philadelphia-data","model_inputs", "phil_flood_hist_categories.csv")))
haz_dict = Dict(haz_cat.year .=> haz_cat.category)

out_dir = joinpath(@__DIR__,"data","shap_runs")
abm_data_files = filter(file -> occursin(r"abm_data.*\.h5$",file), readdir(out_dir))
#for file in abm_data_files
filename = "2011_abm_data_142772.h5"
h5file = h5open(joinpath(out_dir,filename), "r")
pop_dat = h5file["pop_data"]
price_dat = h5file["price_data"]
println("Pop Dataset size: ", size(pop_dat))
println(h5file["pop_vars"][:])
println("Price Dataset size: ", size(price_dat))
println(h5file["price_vars"][:])

#Subset to areas with exposed properties 
flood_year = read(h5file["historical flood year"])
exp_bgs = phil_flood_record[phil_flood_record[!,string(flood_year)] .> 0,"GEOID"]
#Read in damage df 

flpn_index = filter(x -> x !== nothing, indexin(exp_bgs, h5file["GEOID"][:]))

println("Collecting Recovery Data for Flood Event $(flood_year)")

chunk_size = 15950  # Adjust based on your RAM
total_rows = Int(size(pop_dat, 1))

        
println("Processing $total_rows rows in chunks of $chunk_size")


# calculate population trajectories
println("Starting Population Trajectories...")
flpn_low = process_pop(pop_dat;var_col=4, subset=true)
flpn_med = process_pop(pop_dat;var_col=5, subset=true)
flpn_high = process_pop(pop_dat;var_col=6, subset=true)
#Repeat for entire floodplain pop 
pop_flpn = flpn_low .+ flpn_med .+ flpn_high

flpn_norm_low = normal_set(flpn_low)
flpn_norm_med = normal_set(flpn_med)
flpn_norm_high = normal_set(flpn_high)
flpn_norm = normal_set(pop_flpn)

#Save data 
println("Saving Population Data...")
pop_dir = joinpath(out_dir, "post_process","population",haz_dict[flood_year])
mkpath(pop_dir)
low_df = DataFrame(flpn_norm_low,Symbol.(1979 .+ collect(1:size(flpn_norm_low,2))))
CSV.write(joinpath(pop_dir,"$(flood_year)_model_outcome_flpn_pop_norm_low.csv"), low_df)
med_df = DataFrame(flpn_norm_med,Symbol.(1979 .+ collect(1:size(flpn_norm_med,2))))
CSV.write(joinpath(pop_dir,"$(flood_year)_model_outcome_flpn_pop_norm_med.csv"), med_df)
high_df = DataFrame(flpn_norm_high,Symbol.(1979 .+ collect(1:size(flpn_norm_high,2))))
CSV.write(joinpath(pop_dir,"$(flood_year)_model_outcome_flpn_pop_norm_high.csv"), high_df)
pop_df = DataFrame(flpn_norm,Symbol.(1979 .+ collect(1:size(flpn_norm,2))))
CSV.write(joinpath(pop_dir,"$(flood_year)_model_outcome_flpn_pop_norm.csv"), pop_df)

# calculate price trajectories
println("Starting Price Trajectories...")
flpn_price_low = process_price(price_dat;var_col=1, subset=true)
flpn_price_med = process_price(price_dat;var_col=2, subset=true)
flpn_price_high = process_price(price_dat;var_col=3, subset=true)

flpn_np_low = normal_set(flpn_price_low)
flpn_np_med = normal_set(flpn_price_med)
flpn_np_high = normal_set(flpn_price_high)

#Save Data
println("Saving Price Data...")
price_dir = joinpath(out_dir, "post_process","price",haz_dict[flood_year])
mkpath(price_dir)
low_np_df = DataFrame(flpn_np_low,Symbol.(1979 .+ collect(1:size(flpn_np_low,2))))
CSV.write(joinpath(price_dir,"$(flood_year)_model_outcome_flpn_avg_price_norm_low.csv"), low_np_df)
med_np_df = DataFrame(flpn_np_med,Symbol.(1979 .+ collect(1:size(flpn_np_med,2))))
CSV.write(joinpath(price_dir,"$(flood_year)_model_outcome_flpn_avg_price_norm_med.csv"), med_np_df)
high_np_df = DataFrame(flpn_np_high,Symbol.(1979 .+ collect(1:size(flpn_np_high,2))))
CSV.write(joinpath(price_dir,"$(flood_year)_model_outcome_flpn_avg_price_norm_high.csv"), high_np_df)

close(h5file)

println("Finished with $(flood_year)!")
#end




