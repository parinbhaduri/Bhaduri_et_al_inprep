#activate project environment
import Pkg
Pkg.activate(".")
Pkg.instantiate()

### PARALLEL ENSEMBLE RUN ###
include(joinpath(dirname(@__DIR__), "src", "config_parallel.jl"))

### For housing development
@everywhere function phil_build(;flood_rec = synth_flood_record, build_perc=build_perc, no_of_years=no_of_years, start_year=start_year, seed=seed)
    model = phil_model(;flood_rec = flood_rec, no_of_years=Int(no_of_years), start_year=Int(start_year), build_inc_perc=build_perc, seed=seed)      
    return model
end

build_params = Dict(
    :build_perc=>collect(range(0.0,0.5,step=0.05)),
    :no_of_years=>39,
    :start_year=>1981, 
    :seed=>1500
)


adf_build,mdf_build = paramscan(build_params, phil_build; parallel=true, showprogress=true, adata=simul_adata, mdata=simul_mdata, n=39)

### For house pricing
@everywhere function phil_price(;flood_rec = synth_flood_record, price_perc=price_perc, no_of_years=no_of_years, start_year=start_year, seed=seed)
    model = phil_model(;flood_rec = flood_rec, no_of_years=Int(no_of_years), start_year=Int(start_year), price_inc_perc=price_perc, seed=seed)      
    return model
end

price_params = Dict(
    :price_perc=>collect(range(0.0,0.5,step=0.05)),
    :no_of_years=>39,
    :start_year=>1981, 
    :seed=>1500
)

adf_price,mdf_price = paramscan(price_params, phil_price; parallel=true, showprogress=true, adata=simul_adata, mdata=simul_mdata, n=39)


using Plots
using ColorSchemes
include("sim_functions.jl")

build_plots = simul_plot(adf_build, :build_perc; leg = :outertopright)
plot(build_plots..., layout=(3, 2), size = (1100, 1000), plot_title = "Pop. Dynamics when changing Build %")
build_market_plots = simul_market(adf_build,mdf_build, :build_perc; leg = :outertopright)
plot(build_market_plots..., layout=(3, 2), size = (1100, 1000), plot_title = "Market Dynamics when changing Build %")

price_plots = simul_plot(adf_price, :price_perc)
plot(price_plots..., layout=(3, 2), size = (1100, 1000), plot_title = "Pop. Dynamics when changing Price %")
price_market_plots = simul_market(adf_price,mdf_price, :price_perc)
plot(price_market_plots..., layout=(3, 2), size = (1100, 1000), plot_title = "Market Dynamics when changing Price %")