###Simulation study to assess the influence of changes in cofficient value for property area and environmental amenity ###

#activate project environment
import Pkg
Pkg.activate(".")
Pkg.instantiate()



util_coef = Dict(1=> [0, 294707, 130553, 128990, 154887, 72443], 
                 2=> [0, 294707, 130553, 128990, 154887, 72443], 
                 3=> [0, 294707, 130553, 128990, 154887, 72443]
) 

### PARALLEL ENSEMBLE RUN ###
include(joinpath(dirname(@__DIR__), "src", "config_parallel.jl"))

##For Area Utility

area_params_high = Dict(
    :area_l=>300000, :age_l=>130553, :stories_l=>128990, :bath_l=>154887, :env_amen_l=>72443, 
    :area_m=>300000, :age_m=>130553, :stories_m=>128990, :bath_m=>154887, :env_amen_m=>72443, 
    :area_h=>[1e3,1e4,1e5,1e6,1e7, 1e8, 1e9], :age_h=>130553, :stories_h=>128990, :bath_h=>154887, :env_amen_h=>72443,
    :no_of_years=>39,
    :start_year=>1981, 
    :seed=>1500
)

adf_high_area,mdf_high_area = paramscan(area_params_high, phil_util; parallel=true, showprogress=true, adata=simul_adata, mdata = simul_mdata, n=39)

area_params_med = Dict(
    :area_l=>300000, :age_l=>130553, :stories_l=>128990, :bath_l=>154887, :env_amen_l=>72443, 
    :area_m=>[1e3,1e4,1e5,1e6,1e7, 1e8, 1e9], :age_m=>130553, :stories_m=>128990, :bath_m=>154887, :env_amen_m=>72443, 
    :area_h=>300000, :age_h=>130553, :stories_h=>128990, :bath_h=>154887, :env_amen_h=>72443,
    :no_of_years=>39,
    :start_year=>1981, 
    :seed=>1500
)

adf_med_area,mdf_med_area = paramscan(area_params_med, phil_util; parallel=true, showprogress=true, adata=simul_adata, mdata = simul_mdata, n=39)


area_params_low = Dict(
    :area_l=>[1e3,1e4,1e5,1e6,1e7, 1e8, 1e9], :age_l=>130553, :stories_l=>128990, :bath_l=>154887, :env_amen_l=>72443, 
    :area_m=>300000, :age_m=>130553, :stories_m=>128990, :bath_m=>154887, :env_amen_m=>72443, 
    :area_h=>300000, :age_h=>130553, :stories_h=>128990, :bath_h=>154887, :env_amen_h=>72443,
    :no_of_years=>39,
    :start_year=>1981, 
    :seed=>1500
)

adf_low_area,mdf_low_area = paramscan(area_params_low, phil_util; parallel=true, showprogress=true, adata=simul_adata, mdata = simul_mdata, n=39)

## For Env. Amenities ##
env_params_high = Dict(
    :area_l=>300000, :age_l=>130553, :stories_l=>128990, :bath_l=>154887, :env_amen_l=>72443, 
    :area_m=>300000, :age_m=>130553, :stories_m=>128990, :bath_m=>154887, :env_amen_m=>72443, 
    :area_h=>300000, :age_h=>130553, :stories_h=>128990, :bath_h=>154887, :env_amen_h=>[1e3,1e4,1e5,1e6,1e7, 1e8,1e9],
    :no_of_years=>39,
    :start_year=>1981, 
    :seed=>1500
)

adf_high_env,mdf_high_env = paramscan(env_params_high, phil_util; parallel=true, showprogress=true, adata=simul_adata, mdata = simul_mdata, n=39)

env_params_med = Dict(
    :area_l=>300000, :age_l=>130553, :stories_l=>128990, :bath_l=>154887, :env_amen_l=>72443, 
    :area_m=>300000, :age_m=>130553, :stories_m=>128990, :bath_m=>154887, :env_amen_m=>[1e3,1e4,1e5,1e6,1e7, 1e8, 1e9], 
    :area_h=>300000, :age_h=>130553, :stories_h=>128990, :bath_h=>154887, :env_amen_h=>72443,
    :no_of_years=>39,
    :start_year=>1981, 
    :seed=>1500
)

adf_med_env,mdf_med_env = paramscan(env_params_med, phil_util; parallel=true, showprogress=true, adata=simul_adata, mdata = simul_mdata, n=39)


env_params_low = Dict(
    :area_l=>300000, :age_l=>130553, :stories_l=>128990, :bath_l=>154887, :env_amen_l=>[1e3,1e4,1e5,1e6,1e7, 1e8,1e9], 
    :area_m=>300000, :age_m=>130553, :stories_m=>128990, :bath_m=>154887, :env_amen_m=>72443, 
    :area_h=>300000, :age_h=>130553, :stories_h=>128990, :bath_h=>154887, :env_amen_h=>72443,
    :no_of_years=>39,
    :start_year=>1981, 
    :seed=>1500
)

adf_low_env,mdf_low_env = paramscan(env_params_low, phil_util; parallel=true, showprogress=true, adata=simul_adata, mdata = simul_mdata, n=39)

rmprocs(workers())

##Plot
using Plots
using ColorSchemes
include("sim_functions.jl")

#Plot Area results
pop_area_plots = util_plot(adf_low_area, adf_med_area, adf_high_area, "area"; color = palette(:BrBG_7))
plot(pop_area_plots..., layout=(3, 2), size = (1100, 1000), plot_title = "Pop. Dynamics when changing Util Coef: Area")
mark_area_plots = util_market(adf_low_area,adf_med_area,adf_high_area, mdf_low_area,mdf_med_area,mdf_high_area, "area"; color = palette(:BrBG_7))
plot(mark_area_plots..., layout=(3, 2), size = (1100, 1000), plot_title = "Market Dynamics when changing Util Coef: Area")


##Plot Env Amen results
pop_env_plots = util_plot(adf_low_env, adf_med_env, adf_high_env, "env_amen"; color = palette(:BrBG_7))
plot(pop_env_plots..., layout=(3, 2), size = (1100, 1000), plot_title = "Pop. Dynamics when changing Util Coef: Env. Amenities")
mark_env_plots = util_market(adf_low_env,adf_med_env,adf_high_env, mdf_low_env,mdf_med_env,mdf_high_env, "env_amen"; color = palette(:BrBG_7))
plot(mark_env_plots..., layout=(3, 2), size = (1100, 1000), plot_title = "Market Dynamics when changing Util Coef: Env. Amenities")

