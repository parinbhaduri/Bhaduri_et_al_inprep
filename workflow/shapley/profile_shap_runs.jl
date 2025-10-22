##Time running model and collecting data for shapley
#activate project environment
import Pkg
Pkg.activate(".")
Pkg.instantiate()


using Distributed
addprocs(12, exeflags="--project=$(Base.active_project())")

# instantiate and precompile environment
@everywhere begin
  using Pkg;Pkg.activate(".");
  Pkg.instantiate(); Pkg.precompile()
end


### PARALLEL ENSEMBLE RUN ###
@everywhere begin
    using Dates
    using ProgressMeter
    using CSV, DataFrames
    using Statistics
    using DataStructures
    using Agents
    using CHANCE_C
    using LinearAlgebra
    using HDF5
end

using BenchmarkTools, TimerOutputs

tmr = TimerOutput()

@everywhere include(joinpath(dirname(@__DIR__),"src/data_collect.jl"))
@everywhere include("shap_functions.jl")

#Define model run function with timing checkpoints
function test_run_single(
    params::Tuple,
    output_params::Vector{Symbol},
    initialize;
    n = 1,
    adata=shap_adata,
    shock=false, #Whether we're simulating one large flood shock (true) or multiple small (false)
    kwargs...,
)
    output_params_dict = Dict(output_params .=> params)
    @timeit tmr "Model Initialization" model = initialize(;output_params_dict...)
    #Run
    pop_shares_df = DataFrame() 
    if shock
        @timeit tmr "Model Run (0-4)" df_agent_single_1,_ = run!(model, 4; adata=adata,kwargs...) #Run till year before flood shock
        #Collect shares of agents in every block group
        @timeit tmr "Collect Pop Share" pop_shares_df = shock_pop_shares(model)
        #Continue running till end of time horizon
        @timeit tmr "Model Run (5-39)" df_agent_single_2,_ = run!(model, n-4; adata=adata, init=false, kwargs...)

        df_agent_single = vcat(df_agent_single_1, df_agent_single_2)
    else 
        @timeit tmr "Model Run (No Stopping)" df_agent_single,_ = run!(model, n; adata=adata,kwargs...)
    end
    
    
    #Drop rows in queue
    @timeit tmr "Output Cleaning" begin
        queue_pos = df_agent_single[df_agent_single.agent_type .== Symbol("CHANCE_C.Queue"),:].pos[1:2]
        subset!(df_agent_single, :pos => x -> .!(x .∈ Ref(queue_pos)))
        #Remove missing values
        data_df = combine(groupby(df_agent_single,[:time, :pos]),
            :id => minimum => :bg_id,
            :GEOID .=> (col -> minimum(skipmissing(col))) .=> :GEOID,
            Symbol.(adata[3:end]) .=> (col -> sum(skipmissing(col))) .=> (string.(adata[3:end]) .* "_sum")
        )  
    end
    
    return (data_df, pop_shares_df)
end

#Load calibrated parameter combinations
param_path = joinpath(dirname(@__DIR__),"calibration","data/param_comb_final_mean_thresh_6_ens_250.csv")
calib_combs = DataFrame(CSV.File(param_path))[:,1:14]

#Load flood hazard categories
haz_cat = DataFrame(CSV.File(joinpath(dirname(pwd()), "philadelphia-data","model_inputs", "phil_flood_hist_categories.csv")))

events = combine(groupby(haz_cat, "category")) do group
    # Find indices for min and max flood extents
    min_idx = argmin(group.total_extents)
    max_idx = argmax(group.total_extents)
    
    # For median, sort and find middle index
    sorted_indices = sortperm(group.total_extents)
    median_idx = sorted_indices[div(length(sorted_indices) + 1, 2)]
    
    (
        year_min = group.year[min_idx],
        year_med = group.year[median_idx],
        year_max = group.year[max_idx],
        min_extent = group.total_extents[min_idx],
        median_extent = group.total_extents[median_idx],
        max_extent = group.total_extents[max_idx]
    )
end

flood_years = [2011] #vcat(events.year_min,events.year_med, events.year_max)
one_shock = true
repeat_shocks = false

# Set up additional parameters
add_params = OrderedDict(
    :pop_no=>[0,1,2,3,4,5,6,7,8,9,10],
    :seed=>collect(range(1000,1249))
)

# Extract parameter combinations and names
p_combs = collect((Tuple(row) for row in eachrow(calib_combs)))
output_params = collect(Symbol.(names(calib_combs)))
append!(output_params,collect(keys(add_params)),[:flood_event_year, :flood_repeat])

# Generate all combinations outside of flood shock characteristics
mod_combs = [(c..., p, s, flood_years[1], repeat_shocks) for (c, p, s) in Iterators.product(p_combs,values(add_params)...)];

sim_df, pop_df = test_run_single(mod_combs[11], output_params, PhilPopABM; adata=shap_adata, n=39, shock=one_shock)
show(tmr)
reset_timer!(tmr)






filename = 
pop_share_file = 

save_model_data!(filename, idx, sim_df, n_agents, n_years)
save_pop_share_data!(pop_share_file, idx, pop_df)