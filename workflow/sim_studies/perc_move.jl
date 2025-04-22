#activate project environment
import Pkg
Pkg.activate(".")
Pkg.instantiate()

### PARALLEL ENSEMBLE RUN ###
include(joinpath(dirname(@__DIR__), "src", "config_parallel.jl"))

### For Model Population Growth
@everywhere function phil_pop(;flood_rec = phil_flood_record, perc_growth=perc_growth, no_of_years=no_of_years, start_year=start_year, seed=seed)
    model = phil_model(;flood_rec = flood_rec, no_of_years=Int(no_of_years), start_year=Int(start_year), perc_growth=perc_growth, seed=seed)      
    return model
end

pop_params = Dict(
    :perc_growth=>collect(range(0.0, 0.05, step = 0.01)),
    :no_of_years=>39,
    :start_year=>1981, 
    :seed=>1500
)


adf_pop,mdf_pop = paramscan(pop_params, phil_pop; parallel=true, showprogress=true, adata=simul_adata, mdata=simul_mdata, n=39)

### For Baseline Mode Probability
@everywhere function phil_move(;flood_rec = phil_flood_record, base_move=base_move, no_of_years=no_of_years, start_year=start_year, seed=seed)
    model = phil_model(;flood_rec = flood_rec, no_of_years=Int(no_of_years), start_year=Int(start_year), base_move=base_move, seed=seed)      
    return model
end

move_params = Dict(
    :base_move=>collect(range(0.0, 0.05, step = 0.005)),
    :no_of_years=>39,
    :start_year=>1981, 
    :seed=>1500
)


adf_move,mdf_move = paramscan(move_params, phil_move; parallel=true, showprogress=true, adata=simul_adata, mdata=simul_mdata, n=39)


rmprocs(workers())


using Plots
using ColorSchemes
include("sim_functions.jl")

pop_plots = simul_plot(adf_pop, :perc_growth; leg = :outertopright, color = palette(:Set3_11))
plot(pop_plots..., layout=(3, 2), size = (1100, 1000), plot_title = "Pop. Dynamics when changing Pop. Growth")
pop_mark_plots = simul_market(adf_pop,mdf_pop, :perc_growth; leg = :outertopright, color = palette(:Set3_11))
plot(pop_mark_plots..., layout=(3, 2), size = (1100, 1000), plot_title = "Market Dynamics when changing Pop. Growth")

move_plots = simul_plot(adf_move, :base_move; leg = :outertopright, color = palette(:Set3_11))
plot(move_plots..., layout=(3, 2), size = (1100, 1000), plot_title = "Pop. Dynamics when changing Base Move Prob.")
move_mark_plots = simul_market(adf_move,mdf_move, :base_move; leg = :outertopright, color = palette(:Set3_11))
plot(move_mark_plots..., layout=(3, 2), size = (1100, 1000), plot_title = "Market Dynamics when changing Base Move Prob.")