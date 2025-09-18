###Data
# Calculate Flood matrix and Dict for ABM input
data_location = "philadelphia-data/model_inputs"
bg_file = "phil_flood_bg_2019_nomiss_v1.csv"
flood_file = "phil_flood_hist_year.csv"


##Read in Demographic Data
phil_bg = DataFrame(CSV.File(joinpath(dirname(pwd()), data_location, bg_file)))
pop_files = filter(file -> occursin(r"^philly_cbsa_pop.*\.csv$",file), readdir(joinpath(dirname(pwd()), data_location, "pop_files")))

phil_cbsa_pop = DataFrame[]
for file in pop_files
    df = DataFrame(CSV.File(joinpath(dirname(pwd()), data_location, "pop_files",file)))
    #drop missing values
    dropmissing!(df, :NP)
    #For rows with people and negative income, set income to bottom 10%
    inc_bot_10 = quantile(subset(df, [:NP .=> ByRow(>(0)), :adj_income_2019 .=> ByRow(>(0))]).adj_income_2019, [0.10])[1]
    @. df.adj_income_2019 = ifelse.(df.NP > 0 && df.adj_income_2019 <= 0, inc_bot_10, df.adj_income_2019)
    #Record Population Ensemble member
    df[!,:pop_ens] .= parse(Int, match(r"\d+", file).match)

    push!(phil_cbsa_pop, df)
end

phil_cbsa_pop = vcat(phil_cbsa_pop...)

##Read in flood data
phil_flood = DataFrame(CSV.File(joinpath(dirname(pwd()), data_location, flood_file)))
#transform df to correct format
phil_flood_record = unstack(phil_flood, :GEOID, :year, :perc_flood_extent)
#Extra edits
phil_flood_record[!,"1982"] = zeros(size(phil_flood_record)[1])
select!(phil_flood_record, "GEOID", "1981", "1982", Not(["1982", "2019"]), "2019")


###Functions
function PhilPopABM(;bg_df = phil_bg, phil_df = phil_cbsa_pop, f_df = phil_flood_record,
    perc_growth=0.01, flood_coefficient=0.5, risk_averse=0.5, flood_mem=10, base_move=0.01, build_inc_perc=0.01, price_inc_perc=0.01,
    pop_no=0, penalty=0.5, house_budget_mode="rhea", rhea_coef = 0.7, house_budget_perc=0.33, dist_param = [0.5, 0.25, 0.25],
    prop_l=0.5, env_amen_l=0.5, prop_m=0.5, env_amen_m=0.5, prop_h=0.5, env_amen_h=0.5, seed=seed
)
    #Initialize Pop Distribution
    pop_df = phil_df[phil_df.pop_ens .== pop_no, :]

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


#Function to run model instance and collect data)
function run_single(
    params::Tuple,
    output_params::Vector{Symbol},
    initialize;
    n = 1,
    adata=shap_adata,
    kwargs...,
)
    output_params_dict = Dict(output_params .=> params)
    model = initialize(;output_params_dict...)
    df_agent_single,_ = run!(model, n; adata=adata,kwargs...)
   
    #Drop rows in queue
    queue_pos = df_agent_single[df_agent_single.agent_type .== Symbol("CHANCE_C.Queue"),:].pos[1:2]
    subset!(df_agent_single, :pos => x -> .!(x .∈ Ref(queue_pos)))
    #Remove missing values
    data_df =combine(groupby(df_agent_single,[:time, :pos]),
        :id => minimum => :bg_id,
        Symbol.(adata[2:end]) .=> (col -> sum(skipmissing(col))) .=> (string.(adata[2:end]) .* "_sum")
    )
    return data_df
end

#Write data to hdf5 file  
function save_model_data!(filename, run_idx, df, n_agents=755, n_years=40, n_vars=6)
    h5open(filename, "r+") do file
        # Convert DataFrame to 3D array (agents × years × variables)
        # Assuming df has columns: agent_id, year, var1, var2, var3, var4, var5, var6
        
        # Get variable columns (excluding agent_id and year)
        var_cols = names(df)[4:end]  # Adjust based on your structure
        
        # Convert to 3D array
        data_array = zeros(n_agents, n_years, n_vars)  # agents × years × variables
        
        for row in eachrow(df)
            agent_idx = row.bg_id
            year_idx = row.time + 1 #Time starts at 0
            for (var_idx, var_name) in enumerate(var_cols)
                data_array[agent_idx, year_idx, var_idx] = row[var_name]
            end
        end
        
        # Write to HDF5 (run_idx is the current model run number)
        file["data"][run_idx, :, :, :] = data_array
        
    end
end
#=
comb = combs[5]
adf = run_single(comb, output_params, PhilPopABM; adata=shap_adata, n=39)
save_model_data!(filename, 5, adf, n_agents, n_years, n_vars)

comb = combs[end]

output_params_dict = Dict(output_params .=> params)
model = PhilPopABM(;output_params_dict...)
df_agent_single,_= run!(model, 10; adata=shap_adata)
#Drop rows in queue
queue_pos = df_agent_single[df_agent_single.agent_type .== Symbol("CHANCE_C.Queue"),:].pos[1:2]
subset!(df_agent_single, :pos => x -> .!(x .∈ Ref(queue_pos)))
#Remove missing values
combine(groupby(df_agent_single,[:time, :pos]),
    :id => minimum => :bg_id,
    Symbol.(shap_adata[2:end]) .=> (col -> sum(skipmissing(col))) .=> (string.(shap_adata[2:end]) .* "_sum")
)


adf = run_single(comb, output_params, PhilPopABM; adata=shap_adata, n=39)

save_model_data!(filename, 1, adf, 755, 40, 6)
=#