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
#using Distributed
#addprocs(12, exeflags="--project=$(Base.active_project())")
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


calib_params = Dict(
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
    :seed=>collect(range(1000,6912999))
)


#Calculate number of parameter combinations
combs = Iterators.product(values(calib_params)...)
println("Number of Model Runs: ",length(combs))
flush(stdout)


##Run Models. Collect Data
#Deconstructed model run scheme from Agents.jl
println("Runnning Models. Collecting Data...")
adf_calib, mdf_calib = ModelRuns(calib_params)




#Run Models. Collect Data
#adf_calib_o,mdf_calib_o = paramscan(calib_params, PhilABM; include_constants = true,parallel=true, showprogress=true, adata=calib_adata, mdata=calib_mdata, n=39)


#Save Data as df
#CSV.write(joinpath(@__DIR__,"data/model_run_adf_test.csv"), adf_calib)
#CSV.write(joinpath(@__DIR__,"data/model_run_mdf_test.csv"), mdf_calib)
