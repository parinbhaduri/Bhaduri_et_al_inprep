# load environment and packages
import Pkg
Pkg.activate(".")
Pkg.instantiate()

using Distributed
addprocs(12, exeflags="--project=$(Base.active_project())")

@everywhere begin
    using Random # set seed (since this estimator is Monte Carlo-based) and sample
    using CSV, DataFrames
    using Statistics 
    using MLJ
    using EvoTrees
    using ShapML
    using HDF5
    using ProgressMeter
end

include("shap_functions.jl")
Random.seed!(1)

EvoTreeRegressor = @load EvoTreeRegressor pkg=EvoTrees 

##Load Data ##
#=
function read_model_results(filename)
    h5open(filename, "r") do file
        data = read(file, "data")  # Shape: (60500, 755, 40, 6)
        var_names = read(file, "variable_names")
        parameters = read(file, "parameters")
        
        return data, var_names, parameters
    end
end

# Example: Get results for specific run
function get_run_results(filename, run_idx)
    h5open(filename, "r") do file
        run_data = file["data"][run_idx, :, :, :]  # Shape: (755, 40, 6)
        return run_data
    end
end
=#
###Read in data###
out_dir = joinpath(@__DIR__,"data","shap_DESKTOP")
filename = "abm_data_DESKTOP.h5"

h5file = h5open(joinpath(out_dir,filename), "r")
dataset = h5file["data"]
var_names = h5file["variable_names"]
println("Dataset size: ", size(dataset))
println("Data type: ", eltype(dataset))

#Subset Area to Floodplain (Select Block Groups)
phil_bg = DataFrame(CSV.File(joinpath(dirname(pwd()), "philadelphia-data/model_inputs", "phil_flood_bg_2019_nomiss_v1.csv")))
flpn_bgs = phil_bg[phil_bg.perc_flpn_area .> 0,:]
#Grab BGs impacted from the worst flood on record (2011)
bgs_2011 = phil_flood_record[phil_flood_record[!,"2011"] .> 0,:]
#Select block groups in floodplain that were exposed to the worst flood
flpn_bg_ids = unique(flpn_bgs[flpn_bgs.GEOID .∈ Ref(bgs_2011.GEOID),:].GEOID) #unique(flpn_bgs.GEOID)

flpn_index = indexin(flpn_bg_ids, h5file["GEOID"][:])


chunk_size = 6050  # Adjust based on your RAM
total_rows = Int(size(dataset, 1))
        
println("Processing $total_rows rows in chunks of $chunk_size")
function process_output(dataset;chunk_size=6050,total_rows=60500, var_col = 1, subset=false, index=flpn_index)
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
        flpn_low = nothing
        chunk = nothing
        GC.gc()
    end 
    return output
end



flpn_low = process_output(dataset;var_col=1, subset=true)
flpn_med = process_output(dataset;var_col=2, subset=true)
flpn_high = process_output(dataset;var_col=3, subset=true)

close(h5file)

pop_flpn = flpn_low .+ flpn_med .+ flpn_high
low_prop = flpn_low ./ pop_flpn

out_df = DataFrame(pop_flpn,Symbol.(1979 .+ collect(1:size(output,2))))
CSV.write(joinpath(out_dir, "post_process","model_outcomes","model_outcome_flpn_pop.csv"), out_df)

param_values = DataFrame(CSV.File(joinpath(out_dir,"param_runs_shap.csv")))


# sample random features from the overall set
#idx = rand(1:size(param_values, 1), 24200)
out_df = DataFrame(CSV.File(joinpath(out_dir, "post_process","model_outcomes","model_outcome_flpn_low_prop.csv"))) 
features = copy(param_values)#[idx, :]
targets = copy(out_df)#[idx, :]

#=
#Start with random data for now
out_dir = joinpath(dirname(@__DIR__),"calibration","data","results_135842")
filename = "agents_chunk_1.csv"
data = DataFrame(CSV.File(joinpath(out_dir,filename)))

# Create a parameter combination identifier
param_cols = names(data)[8:end]  # parameter columns
time_col = names(data)[1]          # time column  
data_col = names(data)[2]        # data column


# Group by parameter combinations and create wide format

df_wide = combine(groupby(data, param_cols)) do group_df
    # Create a row with parameter values and time series data
    param_values = NamedTuple(zip(Symbol.(param_cols), first(group_df[!, param_cols])))
    time_series = NamedTuple(zip(Symbol.(group_df[!, time_col]), group_df[!, data_col]))
    merge(param_values, time_series)
end


=#
# define function for parallelized Shapley calculation
@everywhere function predict_var(model, data)
    pred = DataFrame(pop_pred = MLJ.predict(model,data))
    return pred
end

# make regression trees and compute Shapley values
function shapley_reg(yrs, features, targets)
    shap_df = DataFrame(feature_name = names(features))
    @showprogress for yr in yrs
        slr_reg_tree = EvoTreeRegressor(nrounds=200, max_depth=5)
        slr_reg_mach = machine(slr_reg_tree, features, targets[:, Symbol(yr)])
        MLJ.fit!(slr_reg_mach, force=true)

        explain = copy(features)
        reference = copy(features)
        shap_out = ShapML.shap(explain = explain, 
                                reference = reference,
                                model = slr_reg_mach,
                                predict_function = predict_var,
                                sample_size = 100,
                                parallel = :samples, 
                                seed = 1)
        shap_grouped = groupby(shap_out, :feature_name)
        shap_summary = combine(shap_grouped, 
            :shap_effect => (x -> mean(abs.(x))))
        rename!(shap_summary, Dict(:shap_effect_function => Symbol("mean_$(yr)")))
        shap_df = innerjoin(shap_df, shap_summary, on=:feature_name)
    end
    return shap_df
end

yrs = 1980:1:2019
shap_df = shapley_reg(yrs, features, targets)
CSV.write(joinpath(out_dir, "post_process","shapley_indices","shap_indices_flpn_low_prop.csv"), shap_df)