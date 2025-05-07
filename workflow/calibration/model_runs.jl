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
include(joinpath(dirname(@__DIR__), "src", "config_parallel.jl"))
### PARALLEL ENSEMBLE RUN ###
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


calib_params = OrderedDict(
    :risk_averse=>[0.3,0.5,0.7],
    :build_inc_perc=>[0.1,0.25, 0.4],
    :price_inc_perc=>[0.1, 0.2], #could reduce to two
    :rhea_coef=>[0.65, 0.75],
    :base_move=>[0.01,0.02,0.03],
    :prop_l=>[0.25, 0.75], 
    :env_amen_l=>[0.25, 0.75], 
    :prop_m=>[0.25, 0.75],
    :env_amen_m=>[0.25, 0.75], 
    :prop_h=>[0.25, 0.75],
    :env_amen_h=>[0.25, 0.75],
    :penalty=>0.5,#[0.25, 0.75],
    :flood_coefficient=>[0.25, 0.75], 
    :seed=>collect(range(1000,1499)) 
)

#Calculate parameter combinations
combs = Iterators.product(values(calib_params)...)
println("Number of Model Runs: ",length(combs))
flush(stdout)

##Run Models. Collect Data
#Deconstructed model run scheme from Agents.jl
println("Runnning Models. Collecting Data...")
output_params = collect(keys(calib_params))

progress = ProgressMeter.Progress(length(combs); enabled = true)

all_data = ProgressMeter.progress_pmap(combs; progress) do comb 
    run_single(comb, output_params, PhilABM; adata=calib_adata, mdata=calib_mdata, n=39)
end;

println("Writing Data to DF...")

adf_calib = DataFrame()
mdf_calib = DataFrame()
for (df1, df2) in all_data
    append!(adf_calib, df1)
    append!(mdf_calib, df2)
end

#Save Data as df
#CSV.write(joinpath(@__DIR__,"data/model_run_adf_test.csv"), adf_calib)
#CSV.write(joinpath(@__DIR__,"data/model_run_mdf_test.csv"), mdf_calib)

#@time combs = dict_list(calib_params)
#length(combs)

#@time stack(Iterators.product(values(calib_params)...);dims=2)







#Run Models. Collect Data
adf_calib,mdf_calib = paramscan(calib_params, phil_model; parallel=true, showprogress=true, adata=calib_adata, mdata=calib_mdata, n=39)

#Save Data as df
println("Saving Data to CSV...")
CSV.write(joinpath(@__DIR__,"data/model_run_adf.csv"), adf_calib)
CSV.write(joinpath(@__DIR__,"data/model_run_mdf.csv"), mdf_calib)


println("Done!")