#activate project environment
import Pkg
Pkg.activate(".")
Pkg.instantiate()

### PARALLEL ENSEMBLE RUN ###
include(joinpath(dirname(@__DIR__), "src", "config_parallel.jl"))

### For cat penalty

penal_params = Dict(
    :flood_rec => phil_flood_record,
    :risk_averse=>0.7,
    :build_inc_perc=>0.01,
    :perc_growth => 0.01,
    :base_move=>0.01,
    #:penalty=>push!(collect(range(0.0,1000,step=100)), 10000000.0),
    :penalty=>[0,0.3,0.5,0.7,0.9,2,5],
    :no_of_years=>39,
    :start_year=>1981, 
    :seed=>1500
)

adf_penal,mdf_penal = paramscan(penal_params, phil_model; parallel=true, showprogress=true, adata=simul_adata, mdata=simul_mdata, n=39)


### For for flood disamenity
disam_params = Dict(
    :flood_rec => phil_flood_record,
    :risk_averse=>0.7,
    :build_inc_perc=>0.01,
    :perc_growth => 0.01,
    :base_move=>0.01,
    :flood_coefficient=>[0,0.1,0.3,0.5,0.7,0.9,10],
    :no_of_years=>39,
    :start_year=>1981, 
    :seed=>1500
)

adf_disam,mdf_disam = paramscan(disam_params, phil_model; parallel=true, showprogress=true, adata=simul_adata, mdata=simul_mdata, n=39)

rmprocs(workers())


using Plots
using ColorSchemes
include("sim_functions.jl")

penal_plots = simul_plot(adf_penal, :penalty; leg = :outertopright, color = palette(:OrRd_7))
Plots.plot(penal_plots..., layout=(3, 2), size = (1100, 1000), plot_title = "Dynamics when changing Penalty")
mark_penal_plots = simul_market(adf_penal,mdf_penal, :penalty; leg = :outertopright, color = palette(:OrRd_8))
Plots.plot(mark_penal_plots..., layout=(3, 2), size = (1100, 1000), plot_title = "Market Dynamics when changing Penalty")

disam_plots = simul_plot(adf_disam, :flood_coefficient, color = palette(:Blues_7))
Plots.plot(disam_plots..., layout=(3, 2), size = (1100, 1000), plot_title = "Dynamics when changing Flood Coef.")
mark_disam_plots = simul_market(adf_disam,mdf_disam, :flood_coefficient; leg = :outertopright, color = palette(:Blues_7))
Plots.plot(mark_disam_plots..., layout=(3, 2), size = (1100, 1000), plot_title = "Market Dynamics when changing Flood Coef.")
