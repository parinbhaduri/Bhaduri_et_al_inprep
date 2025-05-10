###Simulation study to assess the influence of changes in cofficient value for property area and environmental amenity ###

#activate project environment
import Pkg
Pkg.activate(".")
Pkg.instantiate()


### PARALLEL ENSEMBLE RUN ###
include(joinpath(dirname(@__DIR__), "src", "config_parallel.jl"))

##For Property Utility Coef

prop_params_high = Dict(
    :flood_rec => synth_flood_record,
    :build_inc_perc=>0.5,
    :perc_growth => 0.03,
    :risk_averse=>0.1,
    :base_move=>0.04,
    :prop_l=>0.5, :env_amen_l=>0.5, 
    :prop_m=>0.5, :env_amen_m=>0.5, 
    :prop_h=>[0,0.1,0.3,0.5,0.7,0.9,10], :env_amen_h=>0.5,
    :no_of_years=>39,
    :start_year=>1981, 
    :seed=>1500
)

adf_high_prop,mdf_high_prop = paramscan(prop_params_high, phil_model; parallel=true, showprogress=true, adata=simul_adata, mdata = simul_mdata, n=39)

prop_params_med = Dict(
    :flood_rec => synth_flood_record,
    :risk_averse=>0.1,
    :build_inc_perc=>0.5,
    :perc_growth => 0.03,
    :base_move=>0.04,
    :prop_l=>0.5, :env_amen_l=>0.5, 
    :prop_m=>[0,0.1,0.3,0.5,0.7,0.9,10], :env_amen_m=>0.5, 
    :prop_h=>0.5, :env_amen_h=>0.5,
    :no_of_years=>39,
    :start_year=>1981, 
    :seed=>1500
)

adf_med_prop,mdf_med_prop = paramscan(prop_params_med, phil_model; parallel=true, showprogress=true, adata=simul_adata, mdata = simul_mdata, n=39)


prop_params_low = Dict(
    :flood_rec => synth_flood_record,
    :risk_averse=>0.1,
    :build_inc_perc=>0.5,
    :perc_growth => 0.03,
    :base_move=>0.04,
    :prop_l=>[0,0.1,0.3,0.5,0.7,0.9,10], :env_amen_l=>0.5, 
    :prop_m=>0.5, :env_amen_m=>0.5, 
    :prop_h=>0.5, :env_amen_h=>0.5,
    :no_of_years=>39,
    :start_year=>1981, 
    :seed=>1500
)

adf_low_prop,mdf_low_prop = paramscan(prop_params_low, phil_model; parallel=true, showprogress=true, adata=simul_adata, mdata = simul_mdata, n=39)

## For Env. Amenities ##
env_params_high = Dict(
    :flood_rec => synth_flood_record,
    :risk_averse=>0.1,
    :build_inc_perc=>0.5,
    :perc_growth => 0.03,
    :base_move=>0.04,
    :prop_l=>0.5, :env_amen_l=>0.5, 
    :prop_m=>0.5, :env_amen_m=>0.5, 
    :prop_h=>0.5, :env_amen_h=>[0,0.1,0.3,0.5,0.7,0.9,10],
    :no_of_years=>39,
    :start_year=>1981, 
    :seed=>1500
)

adf_high_env,mdf_high_env = paramscan(env_params_high, phil_model; parallel=true, showprogress=true, adata=simul_adata, mdata = simul_mdata, n=39)

env_params_med = Dict(
    :flood_rec => synth_flood_record,
    :risk_averse=>0.1,
    :build_inc_perc=>0.5,
    :perc_growth => 0.03,
    :base_move=>0.04,
    :prop_l=>0.5, :env_amen_l=>0.5, 
    :prop_m=>0.5, :env_amen_m=>[0,0.1,0.3,0.5,0.7,0.9,10], 
    :prop_h=>0.5, :env_amen_h=>0.5,
    :no_of_years=>39,
    :start_year=>1981, 
    :seed=>1500
)

adf_med_env,mdf_med_env = paramscan(env_params_med, phil_model; parallel=true, showprogress=true, adata=simul_adata, mdata = simul_mdata, n=39)


env_params_low = Dict(
    :flood_rec => synth_flood_record,
    :risk_averse=>0.1,
    :build_inc_perc=>0.5,
    :perc_growth => 0.03,
    :base_move=>0.04,
    :prop_l=>0.5, :env_amen_l=>[0,0.1,0.3,0.5,0.7,0.9,10], 
    :prop_m=>0.5, :env_amen_m=>0.5, 
    :prop_h=>0.5, :env_amen_h=>0.5,
    :no_of_years=>39,
    :start_year=>1981, 
    :seed=>1500
)

adf_low_env,mdf_low_env = paramscan(env_params_low, phil_model; parallel=true, showprogress=true, adata=simul_adata, mdata = simul_mdata, n=39)

rmprocs(workers())

##Plot
using Plots
using ColorSchemes
include("sim_functions.jl")

#Plot prop results
pop_prop_plots = util_plot(adf_low_prop, adf_med_prop, adf_high_prop, "prop"; color = palette(:OrRd_7))
plot(pop_prop_plots..., layout=(3, 2), size = (1100, 1000), plot_title = "Pop. Dynamics when changing Util Coef: Property")
mark_prop_plots = util_market(adf_low_prop,adf_med_prop,adf_high_prop, mdf_low_prop,mdf_med_prop,mdf_high_prop, "prop"; color = palette(:OrRd_7))
plot(mark_prop_plots..., layout=(3, 2), size = (1100, 1000), plot_title = "Market Dynamics when changing Util Coef: Property")


##Plot Env Amen results
pop_env_plots = util_plot(adf_low_env, adf_med_env, adf_high_env, "env_amen"; color = palette(:GnBu_7))
plot(pop_env_plots..., layout=(3, 2), size = (1100, 1000), plot_title = "Pop. Dynamics when changing Util Coef: Env. Amenities")
mark_env_plots = util_market(adf_low_env,adf_med_env,adf_high_env, mdf_low_env,mdf_med_env,mdf_high_env, "env_amen"; color = palette(:GnBu_7))
plot(mark_env_plots..., layout=(3, 2), size = (1100, 1000), plot_title = "Market Dynamics when changing Util Coef: Env. Amenities")



###Look at spread when changing all util params
util_params = Dict(
    :prop_l=>[0.3,0.5,0.7], :env_amen_l=>[0.3,0.5,0.7], 
    :prop_m=>[0.3,0.5,0.7], :env_amen_m=>[0.3,0.5,0.7], 
    :prop_h=>[0.3,0.5,0.7], :env_amen_h=>[0.3,0.5,0.7],
    :no_of_years=>39,
    :start_year=>1981, 
    :seed=>1500
)

adf_util,mdf_util = paramscan(util_params, phil_util; parallel=true, showprogress=true, adata=simul_adata, mdata = simul_mdata, n=39)


#Plot Util results
util_plots = simul_plot(adf_util, :prop_l; leg = false, color = palette(:BrBG_3), lim = (5000,15000))
plot(util_plots..., layout=(3, 2), size = (1100, 1000), plot_title = "Pop. Dynamics when changing R_A")
util_mark_plots = simul_market(adf_util,mdf_util, :prop_l; leg = false, color = palette(:BrBG_3), price_lim =(1e5,7e5))
plot(util_mark_plots..., layout=(3, 2), size = (1100, 1000), plot_title = "Market Dynamics when changing R_A")
