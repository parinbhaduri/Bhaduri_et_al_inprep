### Configuration File ###

#import input data 
include("data_include.jl")



##Create Wrapper of ABM Simulator Function
function PhilSim(bg_df, pop_df, flood_df;no_of_years::Int64, start_year::Int64, perc_growth::Float64, flood_coefficient::Float64, 
    risk_averse::Float64, flood_mem::Int64,  base_move::Float64, build_inc_perc::Float64, price_inc_perc::Float64, 
    penalty::Float64, util_coef::Dict, house_budget_mode::String, house_budget_perc::Float64, seed::Int64
) 
    ### Calculate Flood matrix and Dict for ABM input
    f_dict, f_matrix = CHANCE_C.flood_history(flood_df; no_of_years = no_of_years, start_year = start_year)

    model = Simulator(bg_df, pop_df, f_dict, f_matrix, CHANCE_C.model_step!;
        no_of_years = no_of_years, no_hhs_per_agent = 10, house_budget_mode = house_budget_mode,
        house_choice_mode = "flood_mem_utility", grouped = true, group_col = "adj_income_2019", 
        hh_budget_perc = house_budget_perc, bg_cat = Dict(:col =>"income_cat", :group => [1,2,3]),
        cutoff_dict = OrderedDict(1 => [-60000.00,25000.00], 2 =>[25000.00,75000.00], 3 =>[75000.00, 1e7]),
        simple_anova_coefficients = util_coef, flood_coefficient = flood_coefficient, dist_param = [0.3, 0.4, 0.3],
        penalty = penalty, pop_growth_perc = perc_growth, perc_move = base_move, stock_increase_perc = build_inc_perc, 
        price_increase_perc = price_inc_perc, risk_averse = risk_averse, flood_mem = flood_mem, seed = seed
    )

    return model  
end

###Create function to take different parameter combinations as input:
function phil_model(; flood_rec = phil_flood_record, perc_growth=0.01, flood_coefficient=-50000.0, 
    risk_averse=0.5, flood_mem=10, base_move=0.025, build_inc_perc=0.10, price_inc_perc=0.10, 
    penalty=1000.0, house_budget_mode="rhea", house_budget_perc=0.33,
    area_l=600000, age_l=130553, stories_l=128990, bath_l=154887, env_amen_l=72443, 
    area_m=600000, age_m=130553, stories_m=128990, bath_m=154887, env_amen_m=72443, 
    area_h=600000, age_h=130553, stories_h=128990, bath_h=154887, env_amen_h=72443,
    no_of_years=no_of_years, start_year=start_year, seed=seed
)
util_low = [0, area_l, age_l, stories_l, bath_l, env_amen_l]
util_med = [0, area_m, age_m, stories_m, bath_m, env_amen_m]
util_high = [0, area_h, age_h, stories_h, bath_h, env_amen_h]
util_coef = Dict(1=>util_low, 2=>util_med, 3=>util_high)

model = PhilSim(phil_bg, phil_cbsa_base_pop, flood_rec;no_of_years=Int(no_of_years), start_year=Int(start_year), perc_growth=perc_growth,
                 flood_coefficient=flood_coefficient, risk_averse=risk_averse, flood_mem=flood_mem, base_move=base_move,
                 build_inc_perc=build_inc_perc, price_inc_perc=price_inc_perc, penalty=penalty, util_coef=util_coef, seed=Int(seed),
                 house_budget_mode=house_budget_mode, house_budget_perc=house_budget_perc
)
return model
end
    
