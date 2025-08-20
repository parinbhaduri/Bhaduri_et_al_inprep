#activate project environment
import Pkg
Pkg.activate(".")
Pkg.instantiate()

### PARALLEL ENSEMBLE RUN ###
include(joinpath(dirname(@__DIR__), "src", "config_parallel.jl"))

### For cat penalty
@everywhere function phil_budg(; flood_rec = phil_flood_record, house_budget_mode=house_budget_mode, rhea_coef = rhea_coef, no_of_years=no_of_years, start_year=start_year, seed=seed)
    model = phil_model(;flood_rec = flood_rec, no_of_years=Int(no_of_years), start_year=Int(start_year), house_budget_mode=house_budget_mode, rhea_coef = rhea_coef, seed=seed)      
    return model
end

budg_params = Dict(
    :house_budget_mode=>"rhea",
    :rhea_coef=>collect(range(0.60,0.8,step=0.04)),
    :no_of_years=>39,
    :start_year=>1981, 
    :seed=>1500
)

adf_budg,mdf_budg = paramscan(budg_params, phil_budg; parallel=true, showprogress=true, adata=simul_adata, mdata=simul_mdata, n=39)

using Plots
using ColorSchemes
include("sim_functions.jl")

budg_plots = simul_plot(adf_budg, :rhea_coef; leg = :outertopright, color = (palette(:Oranges, 11)))
plot(budg_plots..., layout=(3, 2), size = (1100, 1000), plot_title = "Pop. Dynamics when changing House Budget")
budg_mark_plots = simul_market(adf_budg,mdf_budget, :rhea_coef; leg = :outertopright,color = (palette(:Oranges, 11)))
plot(budg_mark_plots..., layout=(3, 2), size = (1100, 1000), plot_title = "Market Dynamics when House Budget")