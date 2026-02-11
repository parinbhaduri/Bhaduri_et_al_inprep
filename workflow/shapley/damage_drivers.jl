# load environment and packages
import Pkg
Pkg.activate(".")
Pkg.instantiate()

#using Distributed
#addprocs(4, exeflags="--project=$(Base.active_project())")
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
    using Parquet2
    using CategoricalArrays
end

include(joinpath(dirname(@__DIR__),"src","functions.jl"))
Random.seed!(1)

EvoTreeRegressor = @load EvoTreeRegressor pkg=EvoTrees 


out_dir = joinpath(@__DIR__,"data","shap_runs")

phil_damages = DataFrame(CSV.File(joinpath(dirname(pwd()), "philadelphia-data","flood_hazard", "data","phil_flood_dmg_ens.csv")))
phil_exp = DataFrame(CSV.File(joinpath(dirname(pwd()), "philadelphia-data","model_inputs", "phil_flood_hist_year.csv")))

##Define outcome and hazard variables
haz_size = "Medium"
agent_cats = ["low","high","med"]
#event_years = [1981,1991,2018] #High: [1989,1996,2011], Medium: [1981,1991,2018], Low: [1988,2010,2013]

# define function for parallelized Shapley calculation
@everywhere function predict_var(model, data)
    return DataFrame(pop_pred = MLJ.predict(model,data))
end

# make regression trees and compute Shapley values
function shapley_reg(features, target)
    #shap_df = DataFrame(feature_name = names(features))
    println("Fitting Tree...")
    out_reg_tree = EvoTreeRegressor(nrounds=200, max_depth=5);
    out_reg_mach = machine(out_reg_tree, features, target);
    MLJ.fit!(out_reg_mach, force=true)
    
    explain = copy(features)
    reference = copy(features)
    println("Calculating Shapley Indices...")
    shap_out = ShapML.shap(explain = explain, 
                            reference = reference,
                            model = out_reg_mach,
                            predict_function = predict_var,
                            sample_size = 50, #100
                            parallel = :samples, 
                            seed = 1)
    println("Saving Data...")
    shap_grouped = groupby(shap_out, :feature_name)
    shap_summary = combine(shap_grouped, 
        :shap_effect => (x -> mean(abs.(x))))
    rename!(shap_summary, Dict(:shap_effect_function => Symbol("mean_shap")))
    #shap_df = innerjoin(shap_df, shap_summary, on=:feature_name)
    return shap_summary
end



for agent_cat in agent_cats
    # Load simulated results for each event 
    println("Starting shap analysis for $(agent_cat) income group...")
    println("Loading feature array for each event...")
    ds = Parquet2.Dataset(joinpath(out_dir, "post_process","flood_loss",haz_size))
    filtered_files = findall(file -> occursin(Regex("loss_$(agent_cat)"),file), string.(Parquet2.filelist(ds)))
    all_dfs = mapreduce(vcat, enumerate(filtered_files)) do (i,f)
        append!(ds, f)
        df = ds[i] |> Parquet2.select(:model, :DDF, :burden_value) |> DataFrame
        #Subset to sampled features
        sort!(df, :burden_value)
        n = size(df,1)
        # Select indices that span the entire data range 
        indices = range(1, n, length=159500) |> x -> round.(Int, x)
        sampled_data = df[indices,:]
        #Add Event Year
        m = match(r"(\d{4})", string.(Parquet2.filelist(ds))[f])
        fl_year = parse(Int,m.match)
        sampled_data.year = fl_year
        #=
        event_samples = factor_samples.matrix[factor_samples.matrix.year .== Float64(fl_year),:]
        event_df = innerjoin(df, event_samples, on = [:model, :DDF])

        # Record total flood extent within exposed area for each year
        dmg_bgs = unique(phil_damages[phil_damages[!,"naccs_loss_$(fl_year)"] .> 0.0, :bg_id])
        event_extent = sum(phil_exp[(phil_exp.year .== fl_year) .& (phil_exp.GEOID .∈ Ref(dmg_bgs)), :flood_extent])
        replace!(event_df.year, fl_year => event_extent)
        =#
        return sampled_data
    end

    features = transform!(select(all_dfs, Not(:burden_value)), All() .=> categorical, renamecols=false)
    targets = select(all_dfs, :burden_value)[!,1]

    println("Starting Shapley Index calculation...")

    shap_df = shapley_reg(features, targets)
    CSV.write(joinpath(out_dir, "post_process","shapley_indices","$(haz_size)_fld_shap_indices_flpn_burden_$(agent_cat).csv"), shap_df)
    GC.gc()
end