#activate project environment
import Pkg
Pkg.activate(".")
Pkg.instantiate()

using Distributed, SlurmClusterManager

addprocs(SlurmManager())
@everywhere println("hello from $(myid()):$(gethostname())")

# instantiate and precompile environment
@everywhere begin
  using Pkg;Pkg.activate("."); 
  Pkg.instantiate(); Pkg.precompile()
end

### PARALLEL ENSEMBLE RUN ###
@everywhere begin
    using CSV, DataFrames
    using Statistics
    using DataStructures
    using Agents
    using CHANCE_C
    using LinearAlgebra
end

@everywhere include(joinpath(dirname(@__DIR__),"src/data_collect.jl"))
@everywhere include(joinpath(dirname(@__DIR__),"src/config.jl"))


calib_params = Dict(
    :risk_averse=>0.5,#collect(range(0.1,0.9,step=0.2)),
    :flood_mem=>10,#[5,10,15],
    :build_inc_perc=>0.1,#[0.05,0.1,0.25, 0.4],
    :price_inc_perc=>0.1,#[0.1,0.15,0.2],
    :rhea_coef=>0.7,
    :perc_growth=>0.043,
    :base_move=>0.01,#collect(range(0.0, 0.03, step = 0.005)),
    :prop_l=>0.5,#[0.3,0.5,0.7], 
    :env_amen_l=>0.5,#[0.3,0.5,0.7], 
    :prop_m=>0.5,#[0.3,0.5,0.7],
    :env_amen_m=>0.5,#[0.3,0.5,0.7], 
    :prop_h=>0.5,#[0.3,0.5,0.7],
    :env_amen_h=>0.5,#[0.3,0.5,0.7],
    :penalty=>0.5,
    :flood_coefficient=>0.5,#[0.3,0.5,0.7],
    :no_of_years=>39,
    :start_year=>1981, 
    :seed=>collect(range(1000,1999)) 
)

#Calculate number of parameter combinations
println("Number of Model Runs: ",length(Iterators.product(values(calib_params)...)))
flush(stdout)
"""
count = 0
for iter in Iterators.product(values(calib_params)...)
    if count == 5
        break
    else
        println(iter)
        count += 1
    end
end
"""

#Run Models. Collect Data
adf_calib,mdf_calib = paramscan(calib_params, phil_model; parallel=true, showprogress=true, adata=calib_adata, mdata=calib_mdata, n=39)

#Save Data as df
CSV.write(joinpath(@__DIR__,"data/model_run_adf.csv"), adf_calib)
CSV.write(joinpath(@__DIR__,"data/model_run_mdf.csv"), mdf_calib)
