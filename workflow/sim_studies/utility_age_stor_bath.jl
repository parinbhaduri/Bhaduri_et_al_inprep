###Simulation study to assess the influence of changes in cofficient value for property age, no.stories and no. baths ###

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

##For Age Utility

age_params_high = Dict(
    :area_l=>300000, :age_l=>130553, :stories_l=>128990, :bath_l=>154887, :env_amen_l=>72443, 
    :area_m=>300000, :age_m=>130553, :stories_m=>128990, :bath_m=>154887, :env_amen_m=>72443, 
    :area_h=>300000, :age_h=>[1e3,1e4,1e5,1e6,1e7, 1e8, 1e9], :stories_h=>128990, :bath_h=>154887, :env_amen_h=>72443,
    :no_of_years=>39,
    :start_year=>1981, 
    :seed=>1500
)


adf_high_age,mdf_high_age = paramscan(age_params_high, phil_util; parallel=true, showprogress=true, adata=simul_adata, mdata = simul_mdata, n=39)

age_params_med = Dict(
    :area_l=>300000, :age_l=>130553, :stories_l=>128990, :bath_l=>154887, :env_amen_l=>72443, 
    :area_m=>300000, :age_m=>[1e3,1e4,1e5,1e6,1e7, 1e8, 1e9], :stories_m=>128990, :bath_m=>154887, :env_amen_m=>72443, 
    :area_h=>300000, :age_h=>130553, :stories_h=>128990, :bath_h=>154887, :env_amen_h=>72443,
    :no_of_years=>39,
    :start_year=>1981, 
    :seed=>1500
)

adf_med_age,mdf_med_age = paramscan(age_params_med, phil_util; parallel=true, showprogress=true, adata=simul_adata, mdata = simul_mdata, n=39)


age_params_low = Dict(
    :area_l=>300000, :age_l=>[1e3,1e4,1e5,1e6,1e7, 1e8, 1e9], :stories_l=>128990, :bath_l=>154887, :env_amen_l=>72443, 
    :area_m=>300000, :age_m=>130553, :stories_m=>128990, :bath_m=>154887, :env_amen_m=>72443, 
    :area_h=>300000, :age_h=>130553, :stories_h=>128990, :bath_h=>154887, :env_amen_h=>72443,
    :no_of_years=>39,
    :start_year=>1981, 
    :seed=>1500
)

adf_low_age,mdf_low_age = paramscan(age_params_low, phil_util; parallel=true, showprogress=true, adata=simul_adata, mdata = simul_mdata, n=39)

## For No. Stories ##
stories_params_high = Dict(
    :area_l=>300000, :age_l=>130553, :stories_l=>128990, :bath_l=>154887, :env_amen_l=>72443, 
    :area_m=>300000, :age_m=>130553, :stories_m=>128990, :bath_m=>154887, :env_amen_m=>72443, 
    :area_h=>300000, :age_h=>130553, :stories_h=>[1e3,1e4,1e5,1e6,1e7, 1e8, 1e9], :bath_h=>154887, :env_amen_h=>72443,
    :no_of_years=>39,
    :start_year=>1981, 
    :seed=>1500
)

adf_high_stories,mdf_high_stories = paramscan(stories_params_high, phil_util; parallel=true, showprogress=true, adata=simul_adata, mdata = simul_mdata, n=39)

stories_params_med = Dict(
    :area_l=>300000, :age_l=>130553, :stories_l=>128990, :bath_l=>154887, :env_amen_l=>72443, 
    :area_m=>300000, :age_m=>130553, :stories_m=>[1e3,1e4,1e5,1e6,1e7, 1e8, 1e9], :bath_m=>154887, :env_amen_m=>72443, 
    :area_h=>300000, :age_h=>130553, :stories_h=>128990, :bath_h=>154887, :env_amen_h=>72443,
    :no_of_years=>39,
    :start_year=>1981, 
    :seed=>1500
)

adf_med_stories,mdf_med_stories = paramscan(stories_params_med, phil_util; parallel=true, showprogress=true, adata=simul_adata, mdata = simul_mdata, n=39)


stories_params_low = Dict(
    :area_l=>300000, :age_l=>130553, :stories_l=>[1e3,1e4,1e5,1e6,1e7, 1e8, 1e9], :bath_l=>154887, :env_amen_l=>72443, 
    :area_m=>300000, :age_m=>130553, :stories_m=>128990, :bath_m=>154887, :env_amen_m=>72443, 
    :area_h=>300000, :age_h=>130553, :stories_h=>128990, :bath_h=>154887, :env_amen_h=>72443,
    :no_of_years=>39,
    :start_year=>1981, 
    :seed=>1500
)

adf_low_stories,mdf_low_stories = paramscan(stories_params_low, phil_util; parallel=true, showprogress=true, adata=simul_adata, mdata = simul_mdata, n=39)

## For No. Bathrooms ##
bath_params_high = Dict(
    :area_l=>300000, :age_l=>130553, :stories_l=>128990, :bath_l=>154887, :env_amen_l=>72443, 
    :area_m=>300000, :age_m=>130553, :stories_m=>128990, :bath_m=>154887, :env_amen_m=>72443, 
    :area_h=>300000, :age_h=>130553, :stories_h=>128990, :bath_h=>[1e3,1e4,1e5,1e6,1e7, 1e8, 1e9], :env_amen_h=>72443,
    :no_of_years=>39,
    :start_year=>1981, 
    :seed=>1500
)

adf_high_bath,mdf_high_bath = paramscan(bath_params_high, phil_util; parallel=true, showprogress=true, adata=simul_adata, mdata = simul_mdata, n=39)

bath_params_med = Dict(
    :area_l=>300000, :age_l=>130553, :stories_l=>128990, :bath_l=>154887, :env_amen_l=>72443, 
    :area_m=>300000, :age_m=>130553, :stories_m=>128990, :bath_m=>[1e3,1e4,1e5,1e6,1e7, 1e8, 1e9], :env_amen_m=>72443, 
    :area_h=>300000, :age_h=>130553, :stories_h=>128990, :bath_h=>154887, :env_amen_h=>72443,
    :no_of_years=>39,
    :start_year=>1981, 
    :seed=>1500
)

adf_med_bath,mdf_med_bath = paramscan(bath_params_med, phil_util; parallel=true, showprogress=true, adata=simul_adata, mdata = simul_mdata, n=39)


bath_params_low = Dict(
    :area_l=>300000, :age_l=>130553, :stories_l=>128990, :bath_l=>[1e3,1e4,1e5,1e6,1e7, 1e8, 1e9], :env_amen_l=>72443, 
    :area_m=>300000, :age_m=>130553, :stories_m=>128990, :bath_m=>154887, :env_amen_m=>72443, 
    :area_h=>300000, :age_h=>130553, :stories_h=>128990, :bath_h=>154887, :env_amen_h=>72443,
    :no_of_years=>39,
    :start_year=>1981, 
    :seed=>1500
)

adf_low_bath,mdf_low_bath = paramscan(bath_params_low, phil_util; parallel=true, showprogress=true, adata=simul_adata, mdata = simul_mdata, n=39)

rmprocs(workers())

##Plot
using Plots
using ColorSchemes
include("sim_functions.jl")

#Plot Age results
pop_age_plots = util_plot(adf_low_age, adf_med_age, adf_high_age, "age"; color = palette(:BrBG_7))
plot(pop_age_plots..., layout=(3, 2), size = (1100, 1000), plot_title = "Pop. Dynamics when changing Util Coef: Age")
mark_age_plots = util_market(adf_low_age,adf_med_age,adf_high_age, mdf_low_age,mdf_med_age,mdf_high_age, "age"; color = palette(:BrBG_7))
plot(mark_age_plots..., layout=(3, 2), size = (1100, 1000), plot_title = "Market Dynamics when changing Util Coef: Age")


##Plot No. of Stories results
pop_stories_plots = util_plot(adf_low_stories, adf_med_stories, adf_high_stories, "stories"; color = palette(:BrBG_7))
plot(pop_stories_plots..., layout=(3, 2), size = (1100, 1000), plot_title = "Pop. Dynamics when changing Util Coef: No. of Stories")
mark_stories_plots = util_market(adf_low_stories,adf_med_stories,adf_high_stories, mdf_low_stories,mdf_med_stories,mdf_high_stories, "stories"; color = palette(:BrBG_7))
plot(mark_stories_plots..., layout=(3, 2), size = (1100, 1000), plot_title = "Market Dynamics when changing Util Coef: No. of Stories")


##Plot No. of Baths results
pop_bath_plots = util_plot(adf_low_bath, adf_med_bath, adf_high_bath, "bath"; color = palette(:BrBG_7))
plot(pop_bath_plots..., layout=(3, 2), size = (1100, 1000), plot_title = "Pop. Dynamics when changing Util Coef: Bath")
mark_bath_plots = util_market(adf_low_bath,adf_med_bath,adf_high_bath, mdf_low_bath,mdf_med_bath,mdf_high_bath, "bath"; color = palette(:BrBG_7))
plot(mark_bath_plots..., layout=(3, 2), size = (1100, 1000), plot_title = "Market Dynamics when changing Util Coef: Bath")