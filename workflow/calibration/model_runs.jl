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




calib_params_1 = OrderedDict(
    :risk_averse=>[0.3,0.7], #0.5,
    :build_inc_perc=>[0.1, 0.4], #0.25,
    :price_inc_perc=>[0.1, 0.2], #could reduce to two
    :rhea_coef=>[0.65, 0.75],
    :base_move=>[0.01,0.03], #0.02,
    :prop_l=>[0.25, 0.75],
    :env_amen_l=>[0.25, 0.75],
    :prop_m=>[0.25, 0.75],
    :env_amen_m=>[0.25, 0.75],
    :prop_h=>[0.25, 0.75],
    :env_amen_h=>[0.25, 0.75],
    :penalty=>[0.25, 0.75],
    :flood_coefficient=>[0.25, 0.75],
    :seed=>collect(range(1000,1199))
)


calib_params_2 = OrderedDict(
    :risk_averse=>[0.3,0.7], #0.5,
    :build_inc_perc=>[0.1, 0.4], #0.25,
    :price_inc_perc=>[0.1, 0.2], #could reduce to two
    :rhea_coef=>[0.65, 0.75],
    :base_move=>[0.01,0.03], #0.02,
    :prop_l=>[0.25, 0.75],
    :env_amen_l=>[0.25, 0.75],
    :prop_m=>[0.25, 0.75],
    :env_amen_m=>[0.25, 0.75],
    :prop_h=>[0.25, 0.75],
    :env_amen_h=>[0.25, 0.75],
    :penalty=>[0.25, 0.75],#[0.25, 0.75],
    :flood_coefficient=>[0.25, 0.75],
    :seed=>collect(range(1200,1399))
)


calib_params_3 = OrderedDict(
    :risk_averse=>[0.3,0.7], #0.5,
    :build_inc_perc=>[0.1, 0.4], #0.25,
    :price_inc_perc=>[0.1, 0.2], #could reduce to two
    :rhea_coef=>[0.65, 0.75],
    :base_move=>[0.01,0.03], #0.02,
    :prop_l=>[0.25, 0.75],
    :env_amen_l=>[0.25, 0.75],
    :prop_m=>[0.25, 0.75],
    :env_amen_m=>[0.25, 0.75],
    :prop_h=>[0.25, 0.75],
    :env_amen_h=>[0.25, 0.75],
    :penalty=>[0.25, 0.75],#[0.25, 0.75],
    :flood_coefficient=>[0.25, 0.75],
    :seed=>collect(range(1400,1499))
)


#Calculate parameter combinations
total_runs = length(Iterators.product(values(calib_params_1)...)) + length(Iterators.product(values(calib_params_2)...)) + length(Iterators.product(values(calib_params_3)...))
println("Number of Model Runs: ",total_runs)
flush(stdout)


##Run Models. Collect Data
println("Runnning Part 1 (Seeds 1000 to 1199)...")
adf_calib_1, mdf_calib_1 = ModelRuns(calib_params_1)
println("Saving Part 1 to CSV...")
CSV.write(joinpath(@__DIR__,"data/model_run_adf_1000_1199.csv"), adf_calib_1)
CSV.write(joinpath(@__DIR__,"data/model_run_mdf_1000_1199.csv"), mdf_calib_1)
println("Runnning Part 2 (Seeds 1200 to 1399)...")
adf_calib_2, mdf_calib_2 = ModelRuns(calib_params_2)
println("Saving Part 2 to CSV...")
CSV.write(joinpath(@__DIR__,"data/model_run_adf_1200_1399.csv"), adf_calib_2)
CSV.write(joinpath(@__DIR__,"data/model_run_mdf_1200_1399.csv"), mdf_calib_2)
println("Runnning Part 3 (Seeds 1400 to 1499)...")
adf_calib_3, mdf_calib_3 = ModelRuns(calib_params_3)
println("Saving Part 3 to CSV...")
CSV.write(joinpath(@__DIR__,"data/model_run_adf_1400_1499.csv"), adf_calib_3)
CSV.write(joinpath(@__DIR__,"data/model_run_mdf_1400_1499.csv"), mdf_calib_3)


println("Model Runs Complete!")
"""
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
"""
#Save Data as df
#CSV.write(joinpath(@__DIR__,"data/model_run_adf_test.csv"), adf_calib)
#CSV.write(joinpath(@__DIR__,"data/model_run_mdf_test.csv"), mdf_calib)


#@time combs = dict_list(calib_params)
#length(combs)


#@time stack(Iterators.product(values(calib_params)...);dims=2)