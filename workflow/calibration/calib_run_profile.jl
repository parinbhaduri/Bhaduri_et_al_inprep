#Time and Benchmark model run functions for calibration
#activate project environment
import Pkg
Pkg.activate(".")
Pkg.instantiate()

using Distributed
addprocs(12, exeflags="--project=$(Base.active_project())")
@everywhere begin
    using ProgressMeter
    using CSV, DataFrames
    using Statistics
    using DataStructures
    using Agents
    using CHANCE_C
    using LinearAlgebra
end

@everywhere include(joinpath(dirname(@__DIR__),"src/data_collect.jl"))
@everywhere include("calib_functions.jl")

bench_params = Dict(
    :risk_averse=>0.3,#collect(range(0.1,0.9,step=0.2)),
    :build_inc_perc=>0.25,#[0.05,0.1,0.25, 0.4],
    :price_inc_perc=>0.1,#[0.1,0.15,0.2],
    :rhea_coef=>0.7,
    :base_move=>0.01,#collect(range(0.0, 0.03, step = 0.005)),
    :prop_l=>0.51,#[0.3,0.5,0.7],
    :env_amen_l=>0.52,#[0.3,0.5,0.7],
    :prop_m=>0.53,#[0.3,0.5,0.7],
    :env_amen_m=>0.54,#[0.3,0.5,0.7],
    :prop_h=>0.55,#[0.3,0.5,0.7],
    :env_amen_h=>0.56,#[0.3,0.5,0.7],
    :penalty=>0.57,
    :flood_coefficient=>0.58,#[0.3,0.5,0.7],
    :seed=>collect(range(1000,1099))
)

function ModelRuns(calib_params)
    combs = Iterators.product(values(calib_params)...)
    output_params = collect(keys(calib_params))
    progress = ProgressMeter.Progress(length(combs); enabled = true)


    all_data = ProgressMeter.progress_pmap(combs; progress) do comb
        run_single(comb, output_params, PhilABM; adata=calib_adata, mdata=calib_mdata, n=39)
    end
    return all_data
end

data_df = ModelRuns(bench_params)
varinfo(r"data_df")

function ModelRunsMat(calib_params)
    combs = Iterators.product(values(calib_params)...)
    output_params = collect(keys(calib_params))

    init_agent_dataframe
    init_model_dataframe
    progress = ProgressMeter.Progress(length(combs); enabled = true)
    all_data = ProgressMeter.progress_pmap(combs; progress) do comb
        run_single_mem(comb, output_params, PhilABM; adata=calib_adata, mdata=calib_mdata, n=39)
    end

    a_mat = Matrix()
    m_mat = Matrix()
    for (df1, df2) in all_data
        append!(adf, df1)
        append!(mdf, df2)
    end
    return adf, mdf
    return all_data
end

data_mat = ModelRunsMat(bench_params)
varinfo(r"data_mat")

data_mat[37][2] == Matrix(data_df[37][2])

[Symbol.(dataname.(calib_adata)); o_p]
o_p = collect(keys(bench_params))

reinterpret(reshape(data_mat[:][1]))
data_mat[37][1] == data_mat[37][2]
first.(data_mat)

typeof(data_mat)
reduce(vcat,first.(data_mat))
##Performance Measure 
using BenchmarkTools

b = @benchmarkable ModelRuns($bench_params) seconds=1800 evals=1 samples = 10

v1_1_time = run(b)
