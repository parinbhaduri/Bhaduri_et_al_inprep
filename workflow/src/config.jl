### Configuration File ###

#import input data 
include("data_include.jl")



##Create Wrapper of ABM Simulator Function
function PhilSim(bg_df, pop_df, flood_df;no_of_years::Int64, start_year::Int64, perc_growth::Float64, flood_coefficient::Float64, stay_prob::Float64,
    standardization::Union{Bool, String}, risk_averse::Float64, flood_mem::Int64,  base_move::Float64, build_inc_perc::Float64, price_inc_perc::Float64, 
    penalty::Float64, util_coef::Dict, house_budget_mode::String, rhea_coef::Float64, house_budget_perc::Float64, dist_param::Vector, seed::Int64
) 
    model = Simulator(bg_df, pop_df, flood_df, CHANCE_C.evolve!;
        start_year = start_year, no_of_years = no_of_years, no_hhs_per_agent = 10, house_budget_mode = house_budget_mode,
        house_choice_mode = "flood_ind_utility", grouped = true, group_col = "adj_income_2019", 
        hh_budget_perc = house_budget_perc, rhea_coef = rhea_coef, bg_cat = Dict(:col =>"income_cat", :occ_cat => [1,2,3]),
        cutoff_dict = OrderedDict(1 => [-60000.00,25000.00], 2 =>[25000.00,75000.00], 3 =>[75000.00, 1e7]),
        simple_anova_coefficients = util_coef, flood_coefficient = flood_coefficient, dist_param = dist_param,
        penalty = penalty, pop_growth_perc = perc_growth, perc_move = base_move, stock_increase_perc = build_inc_perc, stay_prob = stay_prob, 
        standardization = standardization, price_increase_perc = price_inc_perc, risk_averse = risk_averse, flood_mem = flood_mem, seed = seed
    )

    return model  
end

###Create function to take different parameter combinations as input:
function phil_model(; flood_rec = phil_flood_record, perc_growth=0.01, flood_coefficient=0.5, 
    risk_averse=0.5, flood_mem=10, base_move=0.01, build_inc_perc=0.010, price_inc_perc=0.010, 
    penalty=0.5, house_budget_mode="rhea", rhea_coef = 0.7, house_budget_perc=0.33, dist_param = [0.3, 0.4, 0.3],
    prop_l=0.5, env_amen_l=0.5, prop_m=0.5, env_amen_m=0.5, 
    prop_h=0.5, env_amen_h=0.5, standardization = "normal", stay_prob = 1.0,
    no_of_years=no_of_years, start_year=start_year, seed=seed
)
    util_low = [prop_l, env_amen_l]
    util_med = [prop_m, env_amen_m]
    util_high = [prop_h, env_amen_h]
    util_coef = Dict(1=>util_low, 2=>util_med, 3=>util_high)

    model = PhilSim(phil_bg, phil_cbsa_base_pop, flood_rec;no_of_years=Int(no_of_years), start_year=Int(start_year), perc_growth=perc_growth,
                 standardization = standardization, flood_coefficient=flood_coefficient, risk_averse=risk_averse, flood_mem=flood_mem, base_move=base_move,
                 build_inc_perc=build_inc_perc, price_inc_perc=price_inc_perc, penalty=penalty, util_coef=util_coef, seed=Int(seed), stay_prob = stay_prob,
                 house_budget_mode=house_budget_mode, rhea_coef = rhea_coef, house_budget_perc=house_budget_perc, dist_param = dist_param
    )
    return model
end

###Create function to handle different utility values as input
function phil_util(;prop_l=prop_l, env_amen_l=env_amen_l, 
    prop_m=prop_m, env_amen_m=env_amen_m, 
    prop_h=prop_h, env_amen_h=env_amen_h, standardization = standardization,
    base_move = base_move, no_of_years=no_of_years, start_year=start_year, seed=seed
)
    util_low = [prop_l, env_amen_l]
    util_med = [prop_m, env_amen_m]
    util_high = [prop_h, env_amen_h]
    util_coef = Dict(1=>util_low, 2=>util_med, 3=>util_high)

    model = PhilSim(phil_bg, phil_cbsa_base_pop, phil_flood_record;no_of_years=Int(no_of_years), start_year=Int(start_year), perc_growth=0.01, flood_coefficient=0.5, 
         risk_averse=0.5, flood_mem=10, base_move=base_move, build_inc_perc=0.10, price_inc_perc=0.10, dist_param = [0.3, 0.4, 0.3], stay_prob,
         standardization=standardization, penalty=0.5, util_coef=util_coef, seed=Int(seed), house_budget_mode="rhea", rhea_coef = 0.7, house_budget_perc=0.33
    )
    return model
end
