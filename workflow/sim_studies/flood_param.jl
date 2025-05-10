#activate project environment
import Pkg
Pkg.activate(".")
Pkg.instantiate()


### PARALLEL ENSEMBLE RUN ###
include(joinpath(dirname(@__DIR__), "src", "config_parallel.jl"))

@everywhere function phil_averse(;flood_rec = phil_flood_record, risk_averse=risk_averse, no_of_years=no_of_years, start_year=start_year, seed=seed)
    model = phil_model(;flood_rec = flood_rec, no_of_years=Int(no_of_years), start_year=Int(start_year), risk_averse=risk_averse, seed=seed)      
    return model
end

averse_params = Dict(
    :risk_averse=>collect(range(0.1,0.9,step=0.2)),
    :no_of_years=>39,
    :start_year=>1981, 
    :seed=>1500
)



adf_averse,mdf_averse = paramscan(averse_params, phil_averse; parallel=true, showprogress=true, adata=simul_adata, mdata=simul_mdata, n=39)

@everywhere function phil_mem(;flood_rec = phil_flood_record, flood_mem=flood_mem, no_of_years=no_of_years, start_year=start_year, seed=seed)
    model = phil_model(;flood_rec = flood_rec, no_of_years=Int(no_of_years), start_year=Int(start_year), flood_mem=flood_mem, seed=seed)             
    return model
end

mem_params = Dict(
    :flood_mem=>collect(range(5,40,step=5)),
    :no_of_years=>39,
    :start_year=>1981, 
    :seed=>1500
)
#collect(range(5,40,step=5))

adf_mem,mdf_mem = paramscan(mem_params, phil_mem; parallel=true, showprogress=true, adata=simul_adata, mdata=simul_mdata, n=39)

rmprocs(workers())


using Plots
using ColorSchemes
include("sim_functions.jl")

averse_plots = simul_plot(adf_averse, :risk_averse; leg = :outertopright, color = palette(:BrBG_6), lim = (7000,15000))
plot(averse_plots..., layout=(3, 2), size = (1100, 1000), plot_title = "Pop. Dynamics when changing R_A")
averse_mark_plots = simul_market(adf_averse,mdf_averse, :risk_averse; leg = :outertopright, color = palette(:BrBG_6), price_lim =(1e5,1e6))
plot(averse_mark_plots..., layout=(3, 2), size = (1100, 1000), plot_title = "Market Dynamics when changing R_A")

mem_plots = simul_plot(adf_mem, :flood_mem, color = palette(:Blues_8), lim = (7000,15000))
plot(mem_plots..., layout=(3, 2), size = (1100, 1000), plot_title = "Pop. Dynamics when changing Flood Mem")
mem_mark_plots = simul_market(adf_mem,mdf_mem, :flood_mem; leg = :outertopright, color = palette(:Blues_8), price_lim =(1e5,1e6))
plot(mem_mark_plots..., layout=(3, 2), size = (1100, 1000), plot_title = "Market Dynamics when changing Flood Mem")
