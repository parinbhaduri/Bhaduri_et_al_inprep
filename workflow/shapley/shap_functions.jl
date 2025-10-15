###Data
# Calculate Flood matrix and Dict for ABM input
data_location = "philadelphia-data/model_inputs"
bg_file = "phil_flood_bg_2019_nomiss_v1.csv"
flood_file = "phil_flood_hist_year.csv"


##Read in Demographic Data
phil_bg = DataFrame(CSV.File(joinpath(dirname(pwd()), data_location, bg_file)))
pop_dir = joinpath(dirname(pwd()), data_location, "pop_files")
#pop_files = filter(file -> occursin(r"^philly_cbsa_pop.*\.csv$",file), readdir(joinpath(dirname(pwd()), data_location, "pop_files")))

#phil_cbsa_pop = DataFrame[]

#=
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
=#
##Read in flood data
phil_flood = DataFrame(CSV.File(joinpath(dirname(pwd()), data_location, flood_file)))
#transform df to correct format
phil_flood_record = unstack(phil_flood, :GEOID, :year, :perc_flood_extent)
#Extra edits
phil_flood_record[!,"1982"] = zeros(size(phil_flood_record)[1])
select!(phil_flood_record, "GEOID", "1981", "1982", Not(["1982", "2019"]), "2019")
#Create synthetic flood record of no floods 
synth_mat = zeros(size(phil_flood_record[:,2:end]))
synth_flood_record = DataFrame(hcat(phil_flood_record.GEOID,synth_mat), names(phil_flood_record))
synth_flood_record[!,:GEOID] = convert.(Int64, synth_flood_record[!,:GEOID])

###Functions
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

function PhilPopABM(;bg_df = phil_bg, flood_event_year=1981, flood_repeat=false, flood_freq=5,
    perc_growth=0.01, flood_coefficient=0.5, risk_averse=0.5, flood_mem=10, base_move=0.01, build_inc_perc=0.01, price_inc_perc=0.01,
    penalty=0.5, house_budget_mode="rhea", rhea_coef = 0.7, house_budget_perc=0.33, dist_param = [0.3, 0.4, 0.3], pop_no = 0,
    prop_l=0.5, env_amen_l=0.5, prop_m=0.5, env_amen_m=0.5, prop_h=0.5, env_amen_h=0.5, seed=1200
)
    #Initialize flood record
    f_df = init_flood(;ref_year=flood_event_year, repeat=flood_repeat, freq=flood_freq)
    
    #Initialize Pop Distribution
    pop_df = load_pop(pop_no) #phil_df[phil_df.pop_ens .== pop_no, :]

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

function shock_pop_shares(model::ABM)
    #Collect counts of HHAgents' group and occ_cat counts
    ag_data = DataFrame(
        bg_id = [model[id].bg_id for id in allids(model) if model[id] isa HHAgent],
        group = [model[id].group for id in allids(model) if model[id] isa HHAgent],
        occ_cat = [model[id].occ_cat for id in allids(model) if model[id] isa HHAgent]
    )
    filter!(:bg_id => x -> x > 0, ag_data) #exclude Queues
    cat_counts = combine(
            groupby(ag_data, [:bg_id, :group, :occ_cat]),
            nrow => :count
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
    #Record BG GEOIDs
    result.GEOID = [model[id].GEOID for id in result.bg_id]
    #Format
    select!(result, :bg_id,:GEOID, Not([:bg_id, :GEOID]))

    return result

end

#Function to run model instance and collect data)
function run_single(
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
        df_agent_single_1,_ = run!(model, 5; adata=adata,kwargs...)
        #Collect shares of agents in every block group
        pop_shares_df = shock_pop_shares(model)
        #Continue running till end of time horizon
        df_agent_single_2,_ = run!(model, n-5; adata=adata, init=false, kwargs...)

        df_agent_single = vcat(df_agent_single_1, df_agent_single_2)
    else 
        df_agent_single,_ = run!(model, n; adata=adata,kwargs...)
    end
    
   
    #Drop rows in queue
    queue_pos = df_agent_single[df_agent_single.agent_type .== Symbol("CHANCE_C.Queue"),:].pos[1:2]
    subset!(df_agent_single, :pos => x -> .!(x .∈ Ref(queue_pos)))
    #Remove missing values
    data_df =combine(groupby(df_agent_single,[:time, :pos]),
        :id => minimum => :bg_id,
        :GEOID .=> (col -> minimum(skipmissing(col))) .=> :GEOID,
        Symbol.(adata[3:end]) .=> (col -> sum(skipmissing(col))) .=> (string.(adata[3:end]) .* "_sum")
    )
    return (data_df, pop_shares_df)
end

#Write data to hdf5 file  
function save_model_data!(filename, run_idx, df, n_agents=755, n_years=40)
    h5open(filename, "r+") do file
        # Convert DataFrame to 3D array (agents × years × variables)
        # Assuming df has columns: agent_id, year, var1, var2, var3, var4, var5, var6
        
        # Get variable columns (excluding agent_id,year,GEOID)
        pop_cols = names(df)[5:7]  # Adjust based on your structure
        house_cols = names(df)[8:end]
        
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