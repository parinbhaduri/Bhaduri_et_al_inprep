#activate project environment
import Pkg
Pkg.activate(".")
Pkg.instantiate()

using CSV, DataFrames
using Statistics

include(joinpath(dirname(@__DIR__), "src", "functions.jl"))

#read in simulated output data
simul_outputs = DataFrame(CSV.File(joinpath(@__DIR__,"data/calib_sim_output.csv")))
#Separate by agent categories
param_cols = ["env_amen_l","price_inc_perc","rhea_coef","base_move",
"prop_l","prop_m","env_amen_m","risk_averse","build_inc_perc",
"env_amen_h","flood_coefficient","penalty","prop_h"]


sim_out_low = select(simul_outputs, r"_low")
sim_out_med = select(simul_outputs, r"_med")
sim_out_high = select(simul_outputs, r"_high")

#Create df to hold param combinations
params_df = unique(select(simul_outputs, param_cols))

#Calculate average variance among ensembles
params_df[!, :var_l] = zeros(size(params_df)[1])
params_df[!, :var_m] = zeros(size(params_df)[1])
params_df[!, :var_h] = zeros(size(params_df)[1])


for (i,group) in enumerate(groupby(simul_outputs, param_cols))
    sim_out_low = select(group, r"_low")
    sim_out_med = select(group, r"_med")
    sim_out_high = select(group, r"_high")

    params_df[i, :var_l] = calc_var(sim_out_low)
    params_df[i, :var_m] = calc_var(sim_out_med)
    params_df[i, :var_h] = calc_var(sim_out_high)
end

#Save intermediate df
CSV.write(joinpath(@__DIR__, "data/param_comb_initial.csv"), params_df)


#Calculate errors for each model ensemble
@chain simul_outputs begin
    @groupby(param_cols)
    @select(r"_low")
    @transform()
end