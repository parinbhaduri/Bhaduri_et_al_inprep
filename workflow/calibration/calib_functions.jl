
###Data
# Calculate Flood matrix and Dict for ABM input
data_location = "philadelphia-data/model_inputs"
bg_file = "phil_flood_bg_2019_v1.csv"
pop_file = "philly_cbsa_pop_0.csv"
flood_file = "phil_flood_hist_year.csv"

##Read in Demographic Data
phil_bg = DataFrame(CSV.File(joinpath(dirname(pwd()), data_location, bg_file)))
phil_cbsa_base_pop = DataFrame(CSV.File(joinpath(dirname(pwd()), data_location, "pop_files", pop_file)))
#drop missing values
dropmissing!(phil_cbsa_base_pop, :NP)
#For rows with people and negative income, set income to bottom 10%
inc_bot_10 = quantile(subset(phil_cbsa_base_pop, [:NP .=> ByRow(>(0)), :adj_income_2019 .=> ByRow(>(0))]).adj_income_2019, [0.10])[1]
@. phil_cbsa_base_pop.adj_income_2019 = ifelse.(phil_cbsa_base_pop.NP > 0 && phil_cbsa_base_pop.adj_income_2019 <= 0, inc_bot_10, phil_cbsa_base_pop.adj_income_2019)

##Read in flood data
phil_flood = DataFrame(CSV.File(joinpath(dirname(pwd()), data_location, flood_file)))
#transform df to correct format
phil_flood_record = unstack(phil_flood, :GEOID, :year, :perc_flood_extent)
#Extra edits
phil_flood_record[!,"1982"] = zeros(size(phil_flood_record)[1])
select!(phil_flood_record, "GEOID", "1981", "1982", Not(["1982", "2019"]), "2019")

phil_dict, phil_matrix = CHANCE_C.flood_history(phil_flood_record; no_of_years = Int(39), start_year = Int(1981))


###Functions 
function PhilABM(;bg_df = phil_bg, pop_df = phil_cbsa_base_pop, f_dict = phil_dict, f_matrix = phil_matrix,
    perc_growth=0.01, flood_coefficient=0.5, risk_averse=0.5, flood_mem=10, base_move=0.01, build_inc_perc=0.10, price_inc_perc=0.10, 
    penalty=0.5, house_budget_mode="rhea", rhea_coef = 0.7, house_budget_perc=0.33, dist_param = [0.3, 0.4, 0.3],
    prop_l=0.5, env_amen_l=0.5, prop_m=0.5, env_amen_m=0.5, prop_h=0.5, env_amen_h=0.5, seed=seed
)
    util_low = [prop_l, env_amen_l]
    util_med = [prop_m, env_amen_m]
    util_high = [prop_h, env_amen_h]
    util_coef = Dict(1=>util_low, 2=>util_med, 3=>util_high)

    model = Simulator(bg_df, pop_df, f_dict, f_matrix, CHANCE_C.model_step!;
        no_of_years = Int(39), no_hhs_per_agent = 10, house_budget_mode = house_budget_mode,
        house_choice_mode = "flood_ind_utility", grouped = true, group_col = "adj_income_2019", 
        hh_budget_perc = house_budget_perc, rhea_coef = rhea_coef, bg_cat = Dict(:col =>"income_cat", :group => [1,2,3]),
        cutoff_dict = OrderedDict(1 => [-60000.00,25000.00], 2 =>[25000.00,75000.00], 3 =>[75000.00, 1e7]),
        simple_anova_coefficients = util_coef, flood_coefficient = flood_coefficient, dist_param = dist_param,
        penalty = penalty, pop_growth_perc = perc_growth, perc_move = base_move, stock_increase_perc = build_inc_perc, 
        price_increase_perc = price_inc_perc, risk_averse = risk_averse, flood_mem = flood_mem, seed = Int(seed)
    )

    return model
end


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

#Deconstructed model run scheme from Agents.jl
function ModelRuns(calib_params)
    combs = Iterators.product(values(calib_params)...)
    output_params = collect(keys(calib_params))
    progress = ProgressMeter.Progress(length(combs); enabled = true)

    all_data = ProgressMeter.progress_pmap(combs; progress) do comb 
        run_single(comb, output_params, PhilABM; adata=calib_adata, mdata=calib_mdata, n=39)
    end

    adf = DataFrame()
    mdf = DataFrame()
    for (df1, df2) in all_data
        append!(adf, df1)
        append!(mdf, df2)
    end
    return adf, mdf
end





# This function is taken from DrWatson:
function dict_list(c::Union{Dict, OrderedDict})
    iterable_fields = filter(k -> typeof(c[k]) <: Vector, keys(c))
    non_iterables = setdiff(keys(c), iterable_fields)

    iterable_dict = Dict(iterable_fields .=> getindex.(Ref(c), iterable_fields))
    non_iterable_dict = Dict(non_iterables .=> getindex.(Ref(c), non_iterables))

    vec(map(Iterators.product(values(iterable_dict)...)) do vals
        dd = Dict(keys(iterable_dict) .=> vals)
        if isempty(non_iterable_dict)
            dd
        elseif isempty(iterable_dict)
            non_iterable_dict
        else
            merge(non_iterable_dict, dd)
        end
    end)
end

