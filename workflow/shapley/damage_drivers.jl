# load environment and packages
import Pkg
Pkg.activate(".")
Pkg.instantiate()

using Distributed
addprocs(4, exeflags="--project=$(Base.active_project())")
using Distributed, SlurmClusterManager

addprocs(SlurmManager())

# instantiate and precompile environment
@everywhere begin
  using Pkg;Pkg.activate(".");
  Pkg.instantiate(); Pkg.precompile()
end

@everywhere begin
    using Random # set seed (since this estimator is Monte Carlo-based) and sample
    using CSV, DataFrames
    using Statistics 
    using MLJ
    using EvoTrees
    using ShapML
    using HDF5
    using ProgressMeter
    using Agents
    using CHANCE_C
    using StatsBase
end

include(joinpath(dirname(@__DIR__),"src","functions.jl"))
Random.seed!(1)

EvoTreeRegressor = @load EvoTreeRegressor pkg=EvoTrees 


out_dir = joinpath(@__DIR__,"data","shap_runs")

##Define outcome and hazard variables
outcome = "flood_loss"
haz_size = "High"
agent_cat = "high"

# Load simulated results for each event 
filtered_files = filter(file -> occursin(Regex("burden_$(agent_cat)\\.csv\$"),file), readdir(joinpath(out_dir,"post_process",outcome,haz_size)))
println("Loading feature array for each event...")
all_dfs = [CSV.read(joinpath(joinpath(out_dir,"post_process",outcome,haz_size), file), DataFrame, 
                    select=[:model, :DDF]) 
           for file in filtered_files]

features = vcat(all_dfs...)
#Select a subset of features
sub_ind = sample(1:Int(size(features,1)/2), 1000, replace=false)
features_explain = vcat([df[sub_ind, :] for df in all_dfs]...)

targets = vcat([CSV.read(joinpath(joinpath(out_dir,"post_process",outcome,haz_size), file), DataFrame, 
                    select=[:value]) for file in filtered_files]...
)[!,1]

haz_cat = DataFrame(CSV.File(joinpath(dirname(pwd()), "philadelphia-data","model_inputs", "phil_flood_hist_categories.csv")))
fld_extents = zeros(size(features,1))
fld_extents_explain = zeros(size(features_explain,1))

println("Storing flood extents for each event...")
# Record total flood extent within exposed area for each year
for (event_idx,year) in enumerate([1989,2011]) #1996,
    fld_extent = haz_cat[haz_cat.year .== year, :total_extents]
    fld_extents[((event_idx - 1) * Int(size(features,1)/2) + 1):(event_idx * Int(size(features,1)/2))] .= fld_extent[1]
    fld_extents_explain[((event_idx - 1) * Int(size(features_explain,1)/2) + 1):(event_idx * Int(size(features_explain,1)/2))] .= fld_extent[1]
end

features.fld_extents = fld_extents
features_explain.fld_extents = fld_extents_explain

#targets = copy(vcat(out_df[:,1:21],out_df_2))#[idx, :]
#targets = vcat([CSV.read(joinpath(out_dir,"post_process",outcome,haz_size,file), DataFrame) for file in filtered_files]...)
# define function for parallelized Shapley calculation
@everywhere function predict_var(model, data)
    pred = DataFrame(pop_pred = MLJ.predict(model,data))
    return pred
end

# make regression trees and compute Shapley values
function shapley_reg(features, target, explain)
    shap_df = DataFrame(feature_name = names(features))
    println("Fitting Tree...")
    out_reg_tree = EvoTreeRegressor(nrounds=200, max_depth=5);
    out_reg_mach = machine(out_reg_tree, features, target);
    MLJ.fit!(out_reg_mach, force=true)
    
    #explain = copy(features)
    reference = copy(features)
    println("Calculating Shapley Indices...")
    shap_out = ShapML.shap(explain = explain, 
                            reference = reference,
                            model = out_reg_mach,
                            predict_function = predict_var,
                            sample_size = 50, #100
                            parallel = :samples, 
                            seed = 1)
    println("Saving Data for Year...")
    shap_grouped = groupby(shap_out, :feature_name)
    shap_summary = combine(shap_grouped, 
        :shap_effect => (x -> mean(abs.(x))))
    rename!(shap_summary, Dict(:shap_effect_function => Symbol("mean_$(yr)")))
    shap_df = innerjoin(shap_df, shap_summary, on=:feature_name)
    return shap_df
end
println("Starting Shapley Index calculation...")

shap_df = shapley_reg(features, targets, features_explain)
CSV.write(joinpath(out_dir, "post_process","shapley_indices","$(haz_size)_fld_shap_indices_flpn_burden_$(agent_cat)_sub.csv"), shap_df)