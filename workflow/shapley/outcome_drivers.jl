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
end

Random.seed!(1)

EvoTreeRegressor = @load EvoTreeRegressor pkg=EvoTrees 

##Load Data ##
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

# sample random features from the overall set
idx = rand(1:size(df_wide, 1), 1000)
features = df_wide[idx, 1:14]
targets = df_wide[idx, 15:end]

# define function for parallelized Shapley calculation
@everywhere function predict_var(model, data)
    pred = DataFrame(slr_pred = MLJ.predict(model,data))
    return pred
end

# make regression trees and compute Shapley values
function shapley_reg(yrs, features, targets)
    shap_df = DataFrame(feature_name = names(features))
    for yr in yrs
        slr_reg_tree = EvoTreeRegressor(nrounds=20, max_depth=5)
        slr_reg_mach = machine(slr_reg_tree, features, targets[:, Symbol(yr)])
        MLJ.fit!(slr_reg_mach, force=true)

        explain = copy(features)
        reference = copy(features)
        shap_out = ShapML.shap(explain = explain, 
                                reference = reference,
                                model = slr_reg_mach,
                                predict_function = predict_var,
                                sample_size = 10,
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

yrs = 0:1:39
shap_df = shapley_reg(yrs, features, targets)
save(joinpath(@__DIR__, "..", "output", "shapley", "shapley_indices_$scenario.csv"), shap_df)