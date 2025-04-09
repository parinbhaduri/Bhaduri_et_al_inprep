#activate project environment
import Pkg
Pkg.activate(".")
Pkg.instantiate()

### PARALLEL ENSEMBLE RUN ###
include(joinpath(dirname(@__DIR__), "src", "config_parallel.jl"))

### For cat penalty
@everywhere function phil_budg(; house_budget_mode=house_budget_mode, house_budget_perc=house_budget_perc, no_of_years=no_of_years, start_year=start_year, seed=seed)
    util_coef = Dict(1=> [0, 600000, 130553, 128990, 154887, 72443], 
                2=> [0, 600000, 130553, 128990, 154887, 72443], 
                3=> [0, 600000, 130553, 128990, 154887, 72443]
    )
    
    model = PhilSim(phil_bg, phil_cbsa_base_pop, synth_flood_record;no_of_years=Int(no_of_years), start_year=Int(start_year), perc_growth=0.01, flood_coefficient=-50000.0, 
             risk_averse=0.5, flood_mem=10, base_move=0.025, build_inc_perc=0.10, price_inc_perc=0.10, 
             penalty=1000.0, util_coef=util_coef, seed=Int(seed), house_budget_mode=house_budget_mode, house_budget_perc=house_budget_perc
    )
             
            
    return model
end

budg_params = Dict(
    :house_budget_mode=>"perc",
    :house_budget_perc=>collect(range(0.1,0.6,step=0.05)),
    :no_of_years=>39,
    :start_year=>1981, 
    :seed=>1500
)

adata = [(hh_low, sum, HH), (hh_med, sum, HH), (hh_high, sum, HH) , (occ_low, sum, BG), (occ_med, sum, BG), (occ_high, sum, BG)]

adf_budg,_ = paramscan(budg_params, phil_budg; parallel=true, showprogress=true, adata, n=39)

using Plots
using ColorSchemes
include("sim_functions.jl")

budg_plots = simul_plot(adf_budg, :house_budget_perc; leg = :outertopright)
plot(budg_plots..., layout=(3, 2), size = (1100, 1000), plot_title = "Dynamics when changing Budget")