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

# Load Damage and Flood Extent estimates
phil_damages = DataFrame(CSV.File(joinpath(dirname(pwd()), "philadelphia-data","flood_hazard", "data","phil_flood_dmg_ens.csv")))
phil_damages.sow_ind .+= 1 #Align with :DDF 
phil_exp = DataFrame(CSV.File(joinpath(dirname(pwd()), "philadelphia-data", "model_inputs", "phil_flood_hist_year.csv")))

#
study_bg_file = h5open(joinpath(out_dir,"post_process","flood_loss","study_bg_data.h5"), "r")
bg_dat = study_bg_file["bg_data"]
study_GEOIDs = study_bg_file["GEOID"][:]

# Load Model parameter values
param_values = DataFrame(CSV.File(joinpath(dirname(out_dir),"shap_DESKTOP","param_runs_shap.csv")))
transform!(param_values, ["pop_no","seed"] .=> categorical, renamecols=false) #Set population dist and seed value cols as categorical
param_values.row_num = 1:nrow(param_values)


#= Test
haz_size = "High"
agent_cat = "total"
ds = Parquet2.Dataset(joinpath(out_dir, "post_process","flood_loss",haz_size))
if agent_cat == "total"
    filtered_files = findall(file -> occursin(Regex("loss.parquet"),file), string.(Parquet2.filelist(ds)))
else
    filtered_files = findall(file -> occursin(Regex("loss_$(agent_cat)"),file), string.(Parquet2.filelist(ds)))
end
    
append!(ds, filtered_files[1])
df = ds[1] |> Parquet2.select(:model, :DDF, :burden_value) |> DataFrame
# Join with param values
df = innerjoin(param_values,df, on = :row_num => :model)
# Rank damage realizations by total value
m = match(r"(\d{4})", string.(Parquet2.filelist(ds))[filtered_files[1]]) #Extract Event Year
fl_year = parse(Int,m.match)
tot_dam = combine(groupby(phil_damages,:sow_ind), String("naccs_loss_$(fl_year)") => sum => :total_damages)
sort!(tot_dam,:total_damages)
tot_dam.DDF_order = collect(1:nrow(tot_dam)) ./ nrow(tot_dam)

df = innerjoin(df,tot_dam, on = :DDF => :sow_ind)

dmg_bgs = unique(phil_damages[phil_damages[!,"naccs_loss_$(fl_year)"] .> 0.0, :bg_id])
bg_ind = filter(!isnothing,indexin(dmg_bgs,study_GEOIDs))
total_pop = zeros(size(bg_dat,1))
for ind in bg_ind
    if agent_cat == "total"
        total_pop .+= dropdims(sum(bg_dat[:,ind,7:9],dims=2), dims = 2)
    else
        total_pop .+= bg_dat[:,ind,cat+6]
    end
end
event_pop = DataFrame(:row_num => 1:length(total_pop), :total_pop => total_pop)
df = innerjoin(df,event_pop, on = :row_num)
tot_pop = combine(groupby(select(df, [:seed,:total_pop]),:seed), :total_pop => mean => :total_pop) #calculate avg population grouped by seed value
sort!(tot_pop,:total_pop)
tot_pop.seed_order = collect(1:nrow(tot_pop)) ./ nrow(tot_pop) #Create normalized ordering for each realization

df = innerjoin(df,select(tot_pop,[:seed,:seed_order]), on = :seed)
select!(df, Not([:row_num, :DDF,:total_damages, :total_pop, :seed]))
 # Calculate mean burden for each seed value 
target_means = combine(groupby(df, :seed), 
                :burden_value => mean => :encoded_seed
)

# Join back to dataframe
df = leftjoin(df, target_means, on=:seed)
disallowmissing!(df, :encoded_seed)
# Use encoded_seed as a Float feature instead
select!(df, Not(:seed))
#Subset to sampled features
sort!(df, :burden_value)
n = size(df,1)
# Select indices that span the entire data range 
indices = range(1, n, length=159500) |> x -> round.(Int, x)
sampled_data = df[indices,:]

describe(fdf)
describe(sampled_data)
=#
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

##Define outcome and hazard variables
haz_sizes = ["High", "Medium", "Low"]
agent_cats = ["low","med","high", "total"] #
outcome = "loss" #"burden" or "loss"
for haz_size in haz_sizes
    println("Starting $(haz_size) Hazard Events...")
    for (cat,agent_cat) in enumerate(agent_cats)
        # Load simulated results for each event 
        println("Starting shap analysis for $(agent_cat) income group...")
        println("Loading feature array for each event...")
        ds = Parquet2.Dataset(joinpath(out_dir, "post_process","flood_loss",haz_size))
        if agent_cat == "total"
            filtered_files = findall(file -> occursin(Regex("loss.parquet"),file), string.(Parquet2.filelist(ds)))
        else
            filtered_files = findall(file -> occursin(Regex("loss_$(agent_cat)"),file), string.(Parquet2.filelist(ds)))
        end
        all_dfs = mapreduce(vcat, enumerate(filtered_files)) do (i,f)
            append!(ds, f)
            df = ds[i] |> Parquet2.select(:model, :DDF, Symbol("$(outcome)_value")) |> DataFrame
            # Join with param values
            df = innerjoin(param_values,df, on = :row_num => :model)
            # Rank damage realizations by total value
            m = match(r"(\d{4})", string.(Parquet2.filelist(ds))[f]) #Extract Event Year
            fl_year = parse(Int,m.match)
            tot_dam = combine(groupby(phil_damages,:sow_ind), String("naccs_loss_$(fl_year)") => sum => :total_damages) #calculate total damages by DDF realization
            sort!(tot_dam,:total_damages)
            tot_dam.DDF_order = collect(1:nrow(tot_dam)) ./ nrow(tot_dam) #Create normalized ordering for each realization

            df = innerjoin(df,tot_dam, on = :DDF => :sow_ind)
            # Rank seed realizations by total pop
            #First, calculate total pop
            dmg_bgs = unique(phil_damages[phil_damages[!,"naccs_loss_$(fl_year)"] .> 0.0, :bg_id])
            bg_ind = filter(!isnothing,indexin(dmg_bgs,study_GEOIDs))
            total_pop = zeros(size(bg_dat,1))
            for ind in bg_ind
                if agent_cat == "total"
                    total_pop .+= dropdims(sum(bg_dat[:,ind,7:9],dims=2), dims = 2)
                else
                    total_pop .+= bg_dat[:,ind,cat+6]
                end
            end
            event_pop = DataFrame(:row_num => 1:length(total_pop), :total_pop => total_pop)
            df = innerjoin(df,event_pop, on = :row_num)
            tot_pop = combine(groupby(select(df, [:seed,:total_pop]),:seed), :total_pop => mean => :total_pop) #calculate avg population grouped by seed value
            sort!(tot_pop,:total_pop)
            tot_pop.seed_order = collect(1:nrow(tot_pop)) ./ nrow(tot_pop) #Create normalized ordering for each realization

            df = innerjoin(df,select(tot_pop,[:seed,:seed_order]), on = :seed)

            #Remove Extra Columns
            select!(df, Not([:row_num, :DDF, :total_damages, :total_pop, :seed]))
            #Subset to sampled features
            sort!(df, Symbol("$(outcome)_value"))
            n = size(df,1)
            # Select indices that span the entire data range 
            indices = range(1, n, length=159500) |> x -> round.(Int, x)
            sampled_data = df[indices,:]

            #= Calculate mean burden for each seed value 
            target_means = combine(groupby(sampled_data, :seed), 
                            :burden_value => mean => :encoded_seed
            )
            # Join back to dataframe
            sampled_data = leftjoin(sampled_data, target_means, on=:seed)
            disallowmissing!(sampled_data, :encoded_seed)
            select!(sampled_data, Not([:seed]))
            =#
            # Record total flood extent within exposed area for each year
            event_extent = sum(phil_exp[(phil_exp.year .== fl_year) .& (phil_exp.GEOID .∈ Ref(dmg_bgs)), :flood_extent])
            sampled_data.fld_extent .= event_extent
            
            return sampled_data
        end

        features = select(all_dfs, Not(Symbol("$(outcome)_value")))
        targets = select(all_dfs, Symbol("$(outcome)_value"))[!,1]

        println("Starting Shapley Index calculation...")
        
        shap_df = shapley_reg(features, targets)
        CSV.write(joinpath(out_dir, "post_process","shapley_indices","$(haz_size)_fld_shap_indices_flpn_$(outcome)_$(agent_cat).csv"), shap_df)
        GC.gc()
    end
end
#=
ds = Parquet2.Dataset(joinpath(out_dir, "post_process","flood_loss",haz_size))
filtered_files = findall(file -> occursin(Regex("loss_low"),file), string.(Parquet2.filelist(ds)))
all_dfs = mapreduce(vcat, enumerate(filtered_files)) do (i,f)
    append!(ds, f)
    df = ds[i] |> Parquet2.select(:model, :DDF, :burden_value) |> DataFrame
    # Join with param values
    jdf = innerjoin(param_values,df, on = :row_num => :model)
    # Rank damage realizations by total value
    m = match(r"(\d{4})", string.(Parquet2.filelist(ds))[f]) #Extract Event Year
    fl_year = parse(Int,m.match)
    tot_dam = combine(groupby(phil_damages,:sow_ind), String("naccs_loss_$(fl_year)") => sum => :total_damages) #calculate total damages by DDF realization
    sort!(tot_dam,:total_damages)
    tot_dam.DDF_order = collect(1:nrow(tot_dam)) ./ nrow(tot_dam) #Create normalized ordering for each realization

    fdf = innerjoin(jdf,tot_dam, on = :DDF => :sow_ind)
    select!(fdf, Not([:row_num, :DDF, :total_damages]))
    #Subset to sampled features
    sort!(fdf, :burden_value)
    n = size(fdf,1)
    # Select indices that span the entire data range 
    indices = range(1, n, length=159500) |> x -> round.(Int, x)
    sampled_data = fdf[indices,:]
       
    # Record total flood extent within exposed area for each year
    dmg_bgs = unique(phil_damages[phil_damages[!,"naccs_loss_$(fl_year)"] .> 0.0, :bg_id])
    event_extent = sum(phil_exp[(phil_exp.year .== fl_year) .& (phil_exp.GEOID .∈ Ref(dmg_bgs)), :flood_extent])
    sampled_data.fld_extent .= event_extent
        
    return sampled_data
end
=#

