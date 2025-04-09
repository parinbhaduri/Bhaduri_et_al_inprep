#activate project environment
import Pkg
Pkg.activate(".")
Pkg.instantiate()

#Load Packages
using CSV, DataFrames
using DataStructures
using Agents
using Statistics,StatsBase,Distributions
using CategoricalArrays

using CHANCE_C

#Include necessary functions from other scripts
include(joinpath(dirname(@__DIR__), "src", "config.jl"))
include(joinpath(dirname(@__DIR__), "src", "data_collect.jl"))


#Initialize model
no_of_years = 39
start_year = 1981 
perc_growth = 0.01 
flood_coef = -500000.0
r_a = 0.5
mem = 10  
base_move =  0.025
build_inc_perc =  0.10
price_inc_perc =  0.10
penalty = 100.0
util_coef = Dict(1=> [0, 294707, 130553, 128990, 154887, 72443], 
                 2=> [0, 294707, 130553, 128990, 154887, 72443], 
                 3=> [0, 294707, 130553, 128990, 154887, 72443]) 
seed = 1500

phil_abm = PhilSim(phil_bg, phil_cbsa_base_pop, synth_flood_record;no_of_years=no_of_years, start_year=start_year, perc_growth=perc_growth, flood_coefficient=flood_coef, 
risk_averse=r_a, flood_mem=mem, base_move=base_move, build_inc_perc=build_inc_perc, price_inc_perc=price_inc_perc, penalty=penalty, util_coef=util_coef, seed=seed)

adata = [(hh_low, sum, HH), (hh_med, sum, HH), (hh_high, sum, HH) , (occ_low, sum, BG), (occ_med, sum, BG), (occ_high, sum, BG)]

adf,_ = run!(phil_abm, no_of_years; adata)






### PARALLEL ENSEMBLE RUN ###
include(joinpath(dirname(@__DIR__), "src", "config_parallel.jl"))

@everywhere function phil_util(;area_l=area_l, age_l=age_l, stories_l=stories_l, bath_l=bath_l, env_amen_l=env_amen_l, 
        area_m=area_m, age_m=age_m, stories_m=stories_m, bath_m=bath_m, env_amen_m=env_amen_m, 
        area_h=area_h, age_h=age_h, stories_h=stories_h, bath_h=bath_h, env_amen_h=env_amen_h,
        no_of_years=no_of_years, start_year=start_year, seed=seed
    )
    util_low = [0, area_l, age_l, stories_l, bath_l, env_amen_l]
    util_med = [0, area_m, age_m, stories_m, bath_m, env_amen_m]
    util_high = [0, area_h, age_h, stories_h, bath_h, env_amen_h]
    util_coef = Dict(1=>util_low, 2=>util_med, 3=>util_high)

    model = PhilSim(phil_bg, phil_cbsa_base_pop, synth_flood_record;no_of_years=Int(no_of_years), start_year=Int(start_year), perc_growth=0.01, flood_coefficient=-50000.0, 
             risk_averse=0.5, flood_mem=10, base_move=0.025, build_inc_perc=0.10, price_inc_perc=0.10, 
             penalty=1000.0, util_coef=util_coef, seed=Int(seed), house_budget_mode="rhea", house_budget_perc=0.33
    )
    return model
end

area_params_high = Dict(
    :area_l=>300000, :age_l=>130553, :stories_l=>128990, :bath_l=>154887, :env_amen_l=>72443, 
    :area_m=>300000, :age_m=>130553, :stories_m=>128990, :bath_m=>154887, :env_amen_m=>72443, 
    :area_h=>collect(range(100000,1000000,step=50000)), :age_h=>130553, :stories_h=>128990, :bath_h=>154887, :env_amen_h=>72443,
    :no_of_years=>39,
    :start_year=>1981, 
    :seed=>1500
)

adata = [(hh_high, sum, HH), (occ_high, sum, BG)]

adf_high,_ = paramscan(area_params_high, phil_util; parallel=true, showprogress=true, adata, n=39)

area_params_med = Dict(
    :area_l=>300000, :age_l=>130553, :stories_l=>128990, :bath_l=>154887, :env_amen_l=>72443, 
    :area_m=>collect(range(100000,1000000,step=50000)), :age_m=>130553, :stories_m=>128990, :bath_m=>154887, :env_amen_m=>72443, 
    :area_h=>300000, :age_h=>130553, :stories_h=>128990, :bath_h=>154887, :env_amen_h=>72443,
    :no_of_years=>39,
    :start_year=>1981, 
    :seed=>1500
)

adata = [(hh_med, sum, HH), (occ_med, sum, BG)]

adf_med,_ = paramscan(area_params_med, phil_util; parallel=true, showprogress=true, adata, n=39)


area_params_low = Dict(
    :area_l=>collect(range(100000,1000000,step=50000)), :age_l=>130553, :stories_l=>128990, :bath_l=>154887, :env_amen_l=>72443, 
    :area_m=>300000, :age_m=>130553, :stories_m=>128990, :bath_m=>154887, :env_amen_m=>72443, 
    :area_h=>300000, :age_h=>130553, :stories_h=>128990, :bath_h=>154887, :env_amen_h=>72443,
    :no_of_years=>39,
    :start_year=>1981, 
    :seed=>1500
)

adata = [(hh_low, sum, HH), (occ_low, sum, BG)]

adf_low,_ = paramscan(area_params_low, phil_util; parallel=true, showprogress=true, adata, n=39)

rmprocs(workers())


##Plot
using Plots
using ColorSchemes
plots = []
push!(plots, plot(adf_low.time, adf_low.sum_hh_low_HH, group = adf_low.area_l, legend=:outertopright, palette = palette(:BrBG_11), linewidth=2, title="Low Income Pop."))
push!(plots, plot(adf_low.time, adf_low.sum_occ_low_BG, group = adf_low.area_l, legend=:outertopright, palette = palette(:BrBG_11), linewidth=2, title="Low Income Residents"))

push!(plots, plot(adf_med.time, adf_med.sum_hh_med_HH, group = adf_med.area_m, legend=:outertopright, palette = palette(:BrBG_11), linewidth=2, title="Medium Income Pop."))
push!(plots, plot(adf_med.time, adf_med.sum_occ_med_BG, group = adf_med.area_m, legend=:outertopright, palette = palette(:BrBG_11), linewidth=2, title="Medium Income Residents"))

push!(plots, plot(adf_high.time, adf_high.sum_hh_high_HH, group = adf_high.area_h, legend=:outertopright, palette = palette(:BrBG_11), linewidth=2, title="High Income Pop."))
push!(plots, plot(adf_high.time, adf_high.sum_occ_high_BG, group = adf_high.area_h, legend=:outertopright, palette = palette(:BrBG_11), linewidth=2, title="High Income Residents"))

plot(plots..., layout=(3, 2), size = (1100, 1000), plot_title = "Dynamics when changing Util Coef: Area")




### For Env. Amenities ###
area_params_high = Dict(
    :area_l=>300000, :age_l=>130553, :stories_l=>128990, :bath_l=>154887, :env_amen_l=>72443, 
    :area_m=>300000, :age_m=>130553, :stories_m=>128990, :bath_m=>154887, :env_amen_m=>72443, 
    :area_h=>300000, :age_h=>130553, :stories_h=>128990, :bath_h=>154887, :env_amen_h=>collect(range(50000,500000,step=5000)),
    :no_of_years=>39,
    :start_year=>1981, 
    :seed=>1500
)

adata = [(hh_high, sum, HH), (occ_high, sum, BG)]

adf_high,_ = paramscan(area_params_high, phil_util; parallel=true, showprogress=true, adata, n=39)

area_params_med = Dict(
    :area_l=>300000, :age_l=>130553, :stories_l=>128990, :bath_l=>154887, :env_amen_l=>72443, 
    :area_m=>300000, :age_m=>130553, :stories_m=>128990, :bath_m=>154887, :env_amen_m=>collect(range(50000,500000,step=5000)), 
    :area_h=>300000, :age_h=>130553, :stories_h=>128990, :bath_h=>154887, :env_amen_h=>72443,
    :no_of_years=>39,
    :start_year=>1981, 
    :seed=>1500
)

adata = [(hh_med, sum, HH), (occ_med, sum, BG)]

adf_med,_ = paramscan(area_params_med, phil_util; parallel=true, showprogress=true, adata, n=39)


area_params_low = Dict(
    :area_l=>300000, :age_l=>130553, :stories_l=>128990, :bath_l=>154887, :env_amen_l=>collect(range(50000,500000,step=5000)), 
    :area_m=>300000, :age_m=>130553, :stories_m=>128990, :bath_m=>154887, :env_amen_m=>72443, 
    :area_h=>300000, :age_h=>130553, :stories_h=>128990, :bath_h=>154887, :env_amen_h=>72443,
    :no_of_years=>39,
    :start_year=>1981, 
    :seed=>1500
)

adata = [(hh_low, sum, HH), (occ_low, sum, BG)]

adf_low,_ = paramscan(area_params_low, phil_util; parallel=true, showprogress=true, adata, n=39)

rmprocs(workers())


##Plot
plots = []
push!(plots, plot(adf_low.time, adf_low.sum_hh_low_HH, group = adf_low.env_amen_l, legend=:outertopright, palette = palette(:BrBG_11), linewidth=2, title="Low Income Pop."))
push!(plots, plot(adf_low.time, adf_low.sum_occ_low_BG, group = adf_low.env_amen_l, legend=:outertopright, palette = palette(:BrBG_11), linewidth=2, title="Low Income Residents"))

push!(plots, plot(adf_med.time, adf_med.sum_hh_med_HH, group = adf_med.env_amen_m, legend=:outertopright, palette = palette(:BrBG_11), linewidth=2, title="Medium Income Pop."))
push!(plots, plot(adf_med.time, adf_med.sum_occ_med_BG, group = adf_med.env_amen_m, legend=:outertopright, palette = palette(:BrBG_11), linewidth=2, title="Medium Income Residents"))

push!(plots, plot(adf_high.time, adf_high.sum_hh_high_HH, group = adf_high.env_amen_h, legend=:outertopright, palette = palette(:BrBG_11), linewidth=2, title="High Income Pop."))
push!(plots, plot(adf_high.time, adf_high.sum_occ_high_BG, group = adf_high.env_amen_h, legend=:outertopright, palette = palette(:BrBG_11), linewidth=2, title="High Income Residents"))

plot(plots..., layout=(3, 2), size = (1100, 1000), plot_title = "Dynamics when changing Util Coef: Env. Amenities")
