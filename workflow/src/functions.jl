### File holds functions needed for workflow analysis ###

## ABM INITIALIZATION & DATA COLLECTION ##
function init_flood(;ref_year=1981, repeat=false, freq=5, f_df=synth_flood_record, ref_df=phil_flood_record)
    df = copy(f_df)
    #Grab reference flood year
    flood_event = ref_df[:,string(ref_year)]
    if repeat
        flood_year = string.(collect(1980+freq:freq:2019))
    else
        flood_year = "1985"
    end
    
    df[!,flood_year] .= flood_event

    return df
    
end

function load_pop(pop_value)
    df = DataFrame(CSV.File(joinpath(pop_dir, "philly_cbsa_pop_$pop_value.csv")))
    #drop missing values
    dropmissing!(df, :NP)
    #For rows with people and negative income, set income to bottom 10%
    inc_bot_10 = quantile(subset(df, [:NP .=> ByRow(>(0)), :adj_income_2019 .=> ByRow(>(0))]).adj_income_2019, [0.10])[1]
    @. df.adj_income_2019 = ifelse.(df.NP > 0 && df.adj_income_2019 <= 0, inc_bot_10, df.adj_income_2019)

    return df
end

function PhilPopABM(;bg_df = phil_bg, flood_event_year=1981, flood_shock=false, flood_freq=5,synth_df=synth_flood_record, rl_df=phil_flood_record,
    perc_growth=0.01, flood_coefficient=0.5, risk_averse=0.5, flood_mem=10, base_move=0.01, build_inc_perc=0.01, price_inc_perc=0.01,
    penalty=0.5, house_budget_mode="rhea", rhea_coef = 0.7, house_budget_perc=0.33, dist_param = [0.3, 0.4, 0.3], pop_no = 0,
    prop_l=0.5, env_amen_l=0.5, prop_m=0.5, env_amen_m=0.5, prop_h=0.5, env_amen_h=0.5, seed=1200,
)
    #Initialize flood record
    #fl_df = init_flood(;ref_year=flood_event_year, repeat=flood_repeat, 
    #                    freq=flood_freq, f_df=synth_df, ref_df=rl_df
    #)
                    
    #Initialize Pop Distribution
    pop_df = load_pop(pop_no) #phil_df[phil_df.pop_ens .== pop_no, :]

    util_low = [prop_l, env_amen_l]
    util_med = [prop_m, env_amen_m]
    util_high = [prop_h, env_amen_h]
    util_coef = Dict(1=>util_low, 2=>util_med, 3=>util_high)


    model = Simulator(bg_df, pop_df, rl_df, CHANCE_C.evolve!;
        start_year = Int(1981), no_of_years = Int(39), no_hhs_per_agent = 10, house_budget_mode = house_budget_mode,
        house_choice_mode = "flood_ind_utility", grouped = true, group_col = "adj_income_2019", stay_prob = 1.0,
        hh_budget_perc = house_budget_perc, rhea_coef = rhea_coef, bg_cat = Dict(:col =>"income_cat", :occ_cat => [1,2,3]),
        cutoff_dict = OrderedDict(1 => [-60000.00,25000.00], 2 =>[25000.00,75000.00], 3 =>[75000.00, 1e7]),
        simple_anova_coefficients = util_coef, flood_coefficient = flood_coefficient, dist_param = dist_param,
        standardization = "normal", penalty = penalty, pop_growth_perc = perc_growth, perc_move = base_move, stock_increase_perc = build_inc_perc,
        price_increase_perc = price_inc_perc, risk_averse = risk_averse, flood_mem = flood_mem, flood_shock = flood_shock,
        flood_event_year = flood_event_year,seed = Int(seed)
    )

    return model
end

##Create Wrapper of ABM Simulator Function
function PhilABM(;bg_df = phil_bg, pop_df = phil_cbsa_base_pop, f_df = phil_flood_record,
    perc_growth=0.01, flood_coefficient=0.5, risk_averse=0.5, flood_mem=10, base_move=0.01, build_inc_perc=0.01, price_inc_perc=0.01,
    penalty=0.5, house_budget_mode="rhea", rhea_coef = 0.7, house_budget_perc=0.33, dist_param = [0.3, 0.4, 0.3],
    prop_l=0.5, env_amen_l=0.5, prop_m=0.5, env_amen_m=0.5, prop_h=0.5, env_amen_h=0.5, seed=seed
)
    util_low = [prop_l, env_amen_l]
    util_med = [prop_m, env_amen_m]
    util_high = [prop_h, env_amen_h]
    util_coef = Dict(1=>util_low, 2=>util_med, 3=>util_high)


    model = Simulator(bg_df, pop_df, f_df, CHANCE_C.evolve!;
        start_year = Int(1981), no_of_years = Int(39), no_hhs_per_agent = 10, house_budget_mode = house_budget_mode,
        house_choice_mode = "flood_ind_utility", grouped = true, group_col = "adj_income_2019", stay_prob = 1.0,
        hh_budget_perc = house_budget_perc, rhea_coef = rhea_coef, bg_cat = Dict(:col =>"income_cat", :occ_cat => [1,2,3]),
        cutoff_dict = OrderedDict(1 => [-60000.00,25000.00], 2 =>[25000.00,75000.00], 3 =>[75000.00, 1e7]),
        simple_anova_coefficients = util_coef, flood_coefficient = flood_coefficient, dist_param = dist_param,
        standardization = "normal", penalty = penalty, pop_growth_perc = perc_growth, perc_move = base_move, stock_increase_perc = build_inc_perc,
        price_increase_perc = price_inc_perc, risk_averse = risk_averse, flood_mem = flood_mem, seed = Int(seed)
    )

    return model
end


function PhilSim(bg_df, pop_df, flood_df;no_of_years::Int64, start_year::Int64, perc_growth::Float64, flood_coefficient::Float64, stay_prob::Float64,
    standardization::Union{Bool, String}, risk_averse::Float64, flood_mem::Int64,  base_move::Float64, build_inc_perc::Float64, price_inc_perc::Float64, 
    penalty::Float64, util_coef::Dict, house_budget_mode::String, rhea_coef::Float64, house_budget_perc::Float64, dist_param::Vector, seed::Int64,
    flood_event_year::Int64,flood_shock::Bool
) 
    model = Simulator(bg_df, pop_df, flood_df, CHANCE_C.evolve!;
        start_year = start_year, no_of_years = no_of_years, no_hhs_per_agent = 10, house_budget_mode = house_budget_mode,
        house_choice_mode = "flood_ind_utility", grouped = true, group_col = "adj_income_2019", 
        hh_budget_perc = house_budget_perc, rhea_coef = rhea_coef, bg_cat = Dict(:col =>"income_cat", :occ_cat => [1,2,3]),
        cutoff_dict = OrderedDict(1 => [-60000.00,25000.00], 2 =>[25000.00,75000.00], 3 =>[75000.00, 1e7]),
        simple_anova_coefficients = util_coef, flood_coefficient = flood_coefficient, dist_param = dist_param, penalty = penalty, pop_growth_perc = perc_growth, 
        flood_event_year=flood_event_year,flood_shock=flood_shock, perc_move = base_move, stock_increase_perc = build_inc_perc, stay_prob = stay_prob, 
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
    no_of_years=no_of_years, start_year=start_year, flood_event_year = 1981, flood_shock=false, seed=seed
)
    util_low = [prop_l, env_amen_l]
    util_med = [prop_m, env_amen_m]
    util_high = [prop_h, env_amen_h]
    util_coef = Dict(1=>util_low, 2=>util_med, 3=>util_high)

    model = PhilSim(phil_bg, phil_cbsa_base_pop, flood_rec;no_of_years=Int(no_of_years), start_year=Int(start_year), perc_growth=perc_growth,
                 standardization = standardization, flood_coefficient=flood_coefficient, risk_averse=risk_averse, flood_mem=flood_mem, base_move=base_move,
                 build_inc_perc=build_inc_perc, price_inc_perc=price_inc_perc, penalty=penalty, util_coef=util_coef, seed=Int(seed), stay_prob = stay_prob,
                 house_budget_mode=house_budget_mode, rhea_coef = rhea_coef, house_budget_perc=house_budget_perc, dist_param = dist_param,
                 flood_event_year=flood_event_year,flood_shock=flood_shock
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





## CALIBRATION ##

#Function to run model instance and collect data
function run_single(
    params::Tuple,
    output_params::Vector{Symbol},
    initialize;
    n = 1,
    kwargs...,
)
    output_params_dict = Dict(output_params .=> params)
    model = initialize(;output_params_dict...)
    df_agent_single, df_model_single = run!(model, n; kwargs...)
   
    insertcols!(df_agent_single, output_params_dict...)
    insertcols!(df_model_single, output_params_dict...)
    return (df_agent_single, df_model_single)
end


function obj_err(f_x, z)
    return (f_x - z)^2    
end

function model_error(outputs::Vector, obs::Vector)
    return sum(obj_err.(outputs, obs))
end

function calc_err(df_group; obs = zeros(4))
    "Calculates model discrepancy error among
    ensemble members"
    #n = nrow(df_group)
    #data_matrix = Matrix(df_group)
    #errors = Vector{Float64}(undef, n)
    #idx = 1
    #Calculate average simulated output
    avg_output = mean(Matrix(df_group), dims=1)
    #Calculate model discrepancy
    mod_disc = model_error(vec(avg_output), vec(obs))
    #@inbounds for i in 1:n
    #    errors[idx] = model_error(Vector(view(data_matrix, i, :)), Vector(obs))
    #    idx += 1
    #end
    return mod_disc
end

function calc_var(df_group)
    "Calculates variance in errors among
    ensemble members"
    n = nrow(df_group)
    data_matrix = Matrix(df_group)
    errors = Vector{Float64}(undef, binomial(n, 2))
    idx = 1
    
    @inbounds for i in 1:n-1
        for j in i+1:n
            errors[idx] = model_error(Vector(view(data_matrix, i, :)), Vector(view(data_matrix, j, :)))
            idx += 1
        end
    end
    
    return var(errors)
end

## SHAPLEY ANALYSIS ##

function shock_pop_shares(model::ABM)
    #Collect counts of HHAgents' group and occ_cat counts
    hh_ids = [id for id in allids(model) if model[id] isa HHAgent]
    n = length(hh_ids)

    bg_id = Vector{Int}(undef, n)  
    group = Vector{Int}(undef, n)
    occ_cat = Vector{Int}(undef, n)
    income = Vector{Float64}(undef, n)
    hh_pop = Vector{Float64}(undef, n) #Number of households per agent
    pop = Vector{Float64}(undef, n) #Number of people per agent

    for (i, id) in enumerate(hh_ids)
        agent = model[id]
        bg_id[i] = agent.bg_id
        group[i] = agent.group
        occ_cat[i] = agent.occ_cat
        income[i] = agent.income
        hh_pop[i] = agent.no_hhs_per_agent
        pop[i] = agent.no_hhs_per_agent * agent.hh_size
    end

    ag_data = DataFrame(; bg_id, group, occ_cat, income, hh_pop, pop)

    filter!(:bg_id => x -> x > 0, ag_data) #exclude Queues
    cat_counts = combine(
            groupby(ag_data, [:bg_id, :group, :occ_cat]),
            nrow => :count,
            :income => sum => :TotalIncome,
            :hh_pop => sum => :TotalHHPop,
            :pop => sum => :TotalPop,
    )

    # Create complete grid of all group-occ_cat combinations to ensure consistent size
    categories = [1,2,3]
    block_groups_list = [id for id in allids(model) if model[id] isa BlockGroup]

    complete_grid = crossjoin(
        DataFrame(bg_id = block_groups_list),
        DataFrame(group = categories),
        DataFrame(occ_cat = categories)
    )

    # Left join to fill in missing combinations with 0
    result = leftjoin(complete_grid, cat_counts, on=[:bg_id, :group, :occ_cat])
    result.count = coalesce.(result.count, 0)
    result.TotalIncome = coalesce.(result.TotalIncome, 0.0)
    result.TotalHHPop = coalesce.(result.TotalHHPop, 0.0)
    result.TotalPop = coalesce.(result.TotalPop, 0.0)
    #Record BG GEOIDs
    result.GEOID = [model[id].GEOID for id in result.bg_id]
    #Format
    select!(result, :bg_id,:GEOID, Not([:bg_id, :GEOID]))

    return result

end

#Function to run model instance and collect data)
function run_shap_single(
    params::Tuple,
    output_params::Vector{Symbol},
    initialize;
    n = 1,
    adata=shap_adata,
    shock=false, #Whether we're simulating one large flood shock (true) or multiple small (false)
    kwargs...,
)
    output_params_dict = Dict(output_params .=> params)
    model = initialize(;output_params_dict...)
    #Run
    pop_shares_df = DataFrame() 
    if shock
        df_agent_single_1,_ = run!(model, 4; adata=adata,kwargs...) #Run till year before flood shock
        #Collect shares of agents in every block group
        pop_shares_df = shock_pop_shares(model)
        #Continue running till end of time horizon
        df_agent_single_2,_ = run!(model, n-4; adata=adata, init=false, kwargs...)

        df_agent_single = vcat(df_agent_single_1, df_agent_single_2)
    else 
        df_agent_single,_ = run!(model, n; adata=adata,kwargs...)
    end
    
    
    #Drop rows in queue
    queue_pos = df_agent_single[df_agent_single.agent_type .== Symbol("CHANCE_C.Queue"),:].pos[1:2]
    subset!(df_agent_single, :pos => x -> .!(x .∈ Ref(queue_pos)))
    #Remove missing values
    data_df = combine(groupby(df_agent_single,[:time, :pos]),
        :id => minimum => :bg_id,
        :GEOID .=> (col -> minimum(skipmissing(col))) .=> :GEOID,
        Symbol.(adata[3:end]) .=> (col -> sum(skipmissing(col))) .=> (string.(adata[3:end]) .* "_sum")
    )
   
    return (data_df, pop_shares_df)
end

#Function to run model instance and collect pop share data ONLY
function run_pop_single(
    params::Tuple,
    output_params::Vector{Symbol},
    initialize;
    n = 1,
    adata=shap_adata,
    shock=false, #Whether we're simulating one large flood shock (true) or multiple small (false)
    kwargs...,
)
    output_params_dict = Dict(output_params .=> params)
    model = initialize(;output_params_dict...)
    #Run
    pop_shares_df = DataFrame() 
    if shock
        step!(model,4) #Run till year before flood shock
        #Collect shares of agents in every block group
        pop_shares_df = shock_pop_shares(model)
    end
    return pop_shares_df
end

#Write data to hdf5 file  
function save_model_data!(filename, run_idx, df, n_agents=755, n_years=40)
    h5open(filename, "r+") do file
        # Convert DataFrame to 3D array (agents × years × variables)
        # Assuming df has columns: agent_id, year, var1, var2, var3, var4, var5, var6
        
        # Get variable columns (excluding agent_id,year,GEOID)
        pop_cols = names(df)[5:10]  # Adjust based on your structure
        house_cols = names(df)[11:end]
        
        # Convert to 3D array
        pop_array = zeros(n_agents, n_years, length(pop_cols))  # agents × years ×  pop variables
        price_array = zeros(n_agents, n_years, length(house_cols))

        for row in eachrow(df)
            agent_idx = row.bg_id
            year_idx = row.time + 1 #Time starts at 0
            for (var_idx, var_name) in enumerate(pop_cols)
                #Subset by data type and add to relevant array
                pop_array[agent_idx, year_idx, var_idx] = row[var_name]
            end
            for (var_idx, var_name) in enumerate(house_cols)
                price_array[agent_idx, year_idx, var_idx] = row[var_name]
            end
        end
        
        # Write to HDF5 (run_idx is the current model run number)
        file["pop_data"][run_idx, :, :, :] = pop_array
        file["price_data"][run_idx, :, :, :] = price_array
        
    end
end

function save_pop_share_data!(filename, run_idx, df)
    h5open(filename, "r+") do file
        # Convert DataFrame to 2D array 

        # Write to HDF5 (run_idx is the current model run number)
        file["pop_share_data"][run_idx, :, :] = Matrix(df[!,2:end])
    end
end