# load environment and packages
import Pkg
Pkg.activate(".")
Pkg.instantiate()

#using Distributed
#addprocs(15, exeflags="--project=$(Base.active_project())")
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
end

include(joinpath(dirname(@__DIR__),"src","functions.jl"))
Random.seed!(1)

EvoTreeRegressor = @load EvoTreeRegressor pkg=EvoTrees 


out_dir = joinpath(@__DIR__,"data","shap_runs")

param_values = DataFrame(CSV.File(joinpath(dirname(out_dir),"shap_DESKTOP","param_runs_shap.csv")))

##Define outcome and hazard variables
outcome = "population"
haz_size = "High"

haz_cat = DataFrame(CSV.File(joinpath(dirname(pwd()), "philadelphia-data","model_inputs", "phil_flood_hist_categories.csv")))
fld_extents = zeros(size(param_values,1)*3)

# Record total flood extent within exposed area for each year
for (event_idx,year) in enumerate([1989,1996,2011])
    fld_extent = haz_cat[haz_cat.year .== year, :total_extents]
    fld_extents[((event_idx - 1) * size(param_values,1) + 1):(event_idx * size(param_values,1))] .= fld_extent[1]
end

# Load simulated results for each event 
filtered_files = filter(file -> occursin(r"norm_high.csv$",file), readdir(joinpath(out_dir,"post_process",outcome,haz_size)))
#out_df = DataFrame(CSV.File(joinpath(out_dir, "post_process","population","2011_model_outcome_flpn_pop_norm_low.csv"))) 
#out_df_2 = DataFrame(CSV.File(joinpath(out_dir, "post_process","population","1996_model_outcome_flpn_pop_norm_low.csv")))

features = copy(vcat(param_values, param_values, param_values))#[idx, :]
features.fld_extents = fld_extents

#targets = copy(vcat(out_df[:,1:21],out_df_2))#[idx, :]
targets = vcat([CSV.read(joinpath(out_dir,"post_process",outcome,haz_size,file), DataFrame) for file in filtered_files]...)
# define function for parallelized Shapley calculation
@everywhere function predict_var(model, data)
    pred = DataFrame(pop_pred = MLJ.predict(model,data))
    return pred
end

# make regression trees and compute Shapley values
function shapley_reg(yrs, features, targets)
    shap_df = DataFrame(feature_name = names(features))
    @showprogress for yr in yrs
        println("Fitting Tree...")
        out_reg_tree = EvoTreeRegressor(nrounds=200, max_depth=5);
        out_reg_mach = machine(out_reg_tree, features, targets[:, Symbol(yr)]);
        MLJ.fit!(out_reg_mach, force=true)

        explain = copy(features)
        reference = copy(features)
        println("Calculating Shapley Indices...")
        shap_out = ShapML.shap(explain = explain, 
                                reference = reference,
                                model = out_reg_mach,
                                predict_function = predict_var,
                                sample_size = 100, #100
                                parallel = :samples, 
                                seed = 1)
         println("Saving Data for Year...")
        shap_grouped = groupby(shap_out, :feature_name)
        shap_summary = combine(shap_grouped, 
            :shap_effect => (x -> mean(abs.(x))))
        rename!(shap_summary, Dict(:shap_effect_function => Symbol("mean_$(yr)")))
        shap_df = innerjoin(shap_df, shap_summary, on=:feature_name)
    end
    return shap_df
end

yrs = 1980:1:2000
shap_df = shapley_reg(yrs, features, targets)
CSV.write(joinpath(out_dir, "post_process","shapley_indices","$(haz_size)_fld_shap_indices_flpn_pop_norm_high.csv"), shap_df)
