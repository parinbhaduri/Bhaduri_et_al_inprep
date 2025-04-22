#activate project environment
import Pkg
Pkg.activate(".")
Pkg.instantiate()

### PARALLEL ENSEMBLE RUN ###
include(joinpath(dirname(@__DIR__), "src", "config_parallel.jl"))

### For cat penalty
@everywhere function phil_penal(;flood_rec = phil_flood_record, penalty=penalty, no_of_years=no_of_years, start_year=start_year, seed=seed)
    model = phil_model(;flood_rec = flood_rec, no_of_years=Int(no_of_years), start_year=Int(start_year), penalty=penalty, seed=seed)      
    return model
end

penal_params = Dict(
    #:penalty=>push!(collect(range(0.0,1000,step=100)), 10000000.0),
    :penalty=>[0.0,1e2,1e3,1e4,1e5,1e6,1e8,1e10],
    :no_of_years=>39,
    :start_year=>1981, 
    :seed=>1500
)

adf_penal,mdf_penal = paramscan(penal_params, phil_penal; parallel=true, showprogress=true, adata=simul_adata, mdata=simul_mdata, n=39)


### For for flood disamenity
@everywhere function phil_disam(;flood_rec = phil_flood_record, flood_coefficient=flood_coefficient, no_of_years=no_of_years, start_year=start_year, seed=seed)
    model = phil_model(;flood_rec = flood_rec, no_of_years=Int(no_of_years), start_year=Int(start_year), flood_coefficient=flood_coefficient, seed=seed)      
    return model
end

disam_params = Dict(
    :flood_coefficient=>[0.0,1e2,1e3,1e4,1e5,1e6,1e8,1e10],
    :no_of_years=>39,
    :start_year=>1981, 
    :seed=>1500
)

adf_disam,mdf_disam = paramscan(disam_params, phil_disam; parallel=true, showprogress=true, adata=simul_adata, mdata=simul_mdata, n=39)

rmprocs(workers())


using Plots
using ColorSchemes
include("sim_functions.jl")

penal_plots = simul_plot(adf_penal, :penalty; leg = :outertopright, color = palette(:OrRd_8))
plot(penal_plots..., layout=(3, 2), size = (1100, 1000), plot_title = "Dynamics when changing Penalty")
mark_penal_plots = simul_market(adf_penal,mdf_penal, :penalty; leg = :outertopright, color = palette(:OrRd_8))
plot(mark_penal_plots..., layout=(3, 2), size = (1100, 1000), plot_title = "Market Dynamics when changing Penalty")

disam_plots = simul_plot(adf_disam, :flood_coefficient, color = palette(:GnBu_8))
plot(disam_plots..., layout=(3, 2), size = (1100, 1000), plot_title = "Dynamics when changing Flood Coef.")
mark_disam_plots = simul_market(adf_disam,mdf_disam, :flood_coefficient; leg = :outertopright, color = palette(:GnBu_8))
plot(mark_disam_plots..., layout=(3, 2), size = (1100, 1000), plot_title = "Market Dynamics when changing Flood Coef.")
