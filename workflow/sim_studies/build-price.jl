#activate project environment
import Pkg
Pkg.activate(".")
Pkg.instantiate()

### PARALLEL ENSEMBLE RUN ###
include(joinpath(dirname(@__DIR__), "src", "config_parallel.jl"))

### For housing development
@everywhere function phil_build(;flood_rec = phil_flood_record, build_perc=build_perc, no_of_years=no_of_years, start_year=start_year, seed=seed)
    model = phil_model(;flood_rec = flood_rec, no_of_years=Int(no_of_years), start_year=Int(start_year), build_inc_perc=build_perc,seed=seed)      
    return model
end

build_params = Dict(
    #:build_perc=>[0.05,0.1,0.25, 0.4],
    :build_perc=>[0.01,0.02,0.1,0.2],#collect(range(0.0,0.05,step=0.005)),
    :no_of_years=>39,
    :start_year=>1981, 
    :seed=>1500
)


adf_build,mdf_build = paramscan(build_params, phil_build; parallel=true, showprogress=true, adata=simul_adata, mdata=simul_mdata, n=39)

### For house pricing
@everywhere function phil_price(;flood_rec = phil_flood_record, build_perc = 0.01, price_perc=price_perc, no_of_years=no_of_years, start_year=start_year, seed=seed)
    model = phil_model(;flood_rec = flood_rec, no_of_years=Int(no_of_years), start_year=Int(start_year), build_inc_perc = build_perc, price_inc_perc=price_perc, seed=seed)      
    return model
end

price_params = Dict(
    :build_perc=>0.01,
    :price_perc=>collect(range(0.0,0.05,step=0.005)),
    :no_of_years=>39,
    :start_year=>1981, 
    :seed=>1500
)

adf_price,mdf_price = paramscan(price_params, phil_price; parallel=true, showprogress=true, adata=simul_adata, mdata=simul_mdata, n=39)


using Plots
using ColorSchemes
include("sim_functions.jl")

build_plots = simul_plot(adf_build, :build_perc; lim = (150000,400000), leg = :outertopright, color = palette(:BrBG_11))
plot(build_plots..., layout=(3, 2), size = (1100, 1000), plot_title = "Pop. Dynamics when changing Build%")
#savefig(joinpath(@__DIR__,"sim_data","figures", "sim_study_build_perc_no_price_pop.png"))

build_market_plots = simul_market(adf_build,mdf_build, :build_perc; leg = :outertopright, color = palette(:BrBG_11))
plot(build_market_plots..., layout=(3, 2), size = (1100, 1000), plot_title = "Market Dynamics when changing Build %")
#savefig(joinpath(@__DIR__,"sim_data","figures", "sim_study_build_perc_no_price_mark.png"))


price_plots = simul_plot(adf_price, :price_perc; lim = (150000,400000), color = palette(:BrBG_11))
plot(price_plots..., layout=(3, 2), size = (1100, 1000), plot_title = "Pop. Dynamics when changing Price %")
price_market_plots = simul_market(adf_price,mdf_price, :price_perc; price_lim =(1e5,2e6), color = palette(:BrBG_11))
plot(price_market_plots..., layout=(3, 2), size = (1100, 1000), plot_title = "Market Dynamics when changing Price %")



## Look at effects of the combination of the two
@everywhere function phil_bp(;flood_rec = synth_flood_record, build_perc=build_perc, price_perc=price_perc, no_of_years=no_of_years, start_year=start_year, seed=seed)
    model = phil_model(;flood_rec = flood_rec, no_of_years=Int(no_of_years), start_year=Int(start_year), build_inc_perc=build_perc, price_inc_perc=price_perc, seed=seed)      
    return model
end

bp_params = Dict(
    :build_perc=>collect(range(0.0,0.5,step=0.2)),
    :price_perc=>collect(range(0.0,0.5,step=0.2)),
    :no_of_years=>39,
    :start_year=>1981, 
    :seed=>1500
)


adf_bp,mdf_bp = paramscan(bp_params, phil_bp; parallel=true, showprogress=true, adata=simul_adata, mdata=simul_mdata, n=39)

bp_plots = simul_plot(adf_bp, :build_perc; leg = :outertopright, color=palette(:BuPu_3))
plot(bp_plots..., layout=(3, 2), size = (1100, 1000), plot_title = "Pop. Dynamics when changing Build % & Price %")
bp_market_plots = simul_market(adf_bp,mdf_bp, :build_perc; leg = :outertopright, color=palette(:BuPu_3))
plot(bp_market_plots..., layout=(3, 2), size = (1100, 1000), plot_title = "Market Dynamics when changing Build % & Price %")

bp_plots = simul_plot(adf_bp, :price_perc; leg = :outertopright, color=palette(:BuPu_3))
plot(bp_plots..., layout=(3, 2), size = (1100, 1000), plot_title = "Pop. Dynamics when changing Build % & Price %")
bp_market_plots = simul_market(adf_bp,mdf_bp, :price_perc; leg = :outertopright, color=palette(:BuPu_3))
plot(bp_market_plots..., layout=(3, 2), size = (1100, 1000), plot_title = "Market Dynamics when changing Build % & Price %")

