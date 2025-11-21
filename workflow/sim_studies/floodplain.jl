#activate project environment
import Pkg
Pkg.activate(".")
Pkg.instantiate()


### PARALLEL ENSEMBLE RUN ###
include(joinpath(dirname(@__DIR__), "src", "config_parallel.jl"))

#Create alternate Flood scenarios
@everywhere begin
    flood_2011 = init_flood(;ref_year=2011, repeat=false, freq=1)
end

@everywhere function phil_flpn(;flood_rec = phil_flood_record, flood_event_year=2011,flood_shock=true, build_perc = 0.01, price_perc=price_perc, no_of_years=no_of_years, start_year=start_year, seed=seed)
    model = phil_model(;flood_rec = flood_rec, no_of_years=Int(no_of_years), flood_event_year=flood_event_year,flood_shock=flood_shock, start_year=Int(start_year), build_inc_perc=build_perc, price_inc_perc=price_perc, seed=seed)      
    return model
end

flpn_params = Dict(
    :build_perc=>0.01,
    :price_perc=>[0.01, 0.015, 0.02],
    :no_of_years=>39,
    :start_year=>1981, 
    :seed=>1200
)

adf_flpn,mdf_flpn = paramscan(flpn_params, phil_flpn; parallel=true, showprogress=true, adata=simul_adata, mdata=simul_mdata, n=39)


rmprocs(workers())


using Plots
using ColorSchemes
include("sim_functions.jl")

flpn_plots = simul_plot(adf_flpn, :price_perc; leg = :outertopright, color = palette(:BrBG_6), lim = (150000,400000))
Plots.plot(flpn_plots..., layout=(3, 2), size = (1100, 1000), plot_title = "Pop. Dynamics when changing price factor")
flpn_mark_plots = simul_market(adf_flpn,mdf_flpn, :price_perc; leg = :outertopright, color = palette(:BrBG_6), lim = (0,1000), price_lim =(1e5,1e6))
Plots.plot(flpn_mark_plots..., layout=(3, 2), size = (1100, 1000), plot_title = "Market Dynamics when changing price factor")



## PRICE Investigation ##
#Lets collect components that go into pricing calculation and recreate the calculation 
param_path = joinpath(dirname(@__DIR__),"calibration","data/param_comb_final_mean_thresh_6_ens_250.csv")
calib_combs = DataFrame(CSV.File(param_path))[:,1:14]
p_combs = collect((Tuple(row) for row in eachrow(calib_combs)))
output_params = collect(Symbol.(names(calib_combs)))
combs = [(c..., 2011, true) for c in p_combs];
append!(output_params, [:flood_event_year, :flood_shock])

input_params_dict = Dict(output_params .=> combs[1])
model = PhilPopABM(;input_params_dict...)
df_price,_ = run!(model, 20; adata=shap_price_adata)
#Drop rows in queue
queue_pos = df_price[df_price.agent_type .== Symbol("CHANCE_C.Queue"),:].pos[1:2]
subset!(df_price, :pos => x -> .!(x .∈ Ref(queue_pos)))
#Calculate total housing value
transform!(df_price,
    [:price_low, :cap_low] => ByRow(*) => :price_low,
    [:price_med, :cap_med] => ByRow(*) => :price_med,
    [:price_high, :cap_high] => ByRow(*) => :price_high,
)
#Remove missing values
price_df = combine(groupby(df_price,[:time, :pos]),
    :id => minimum => :bg_id,
    :GEOID .=> (col -> minimum(skipmissing(col))) .=> :GEOID,
    Symbol.(shap_price_adata[3:end]) .=> (col -> sum(skipmissing(col))) .=> (string.(shap_price_adata[3:end]))
)

#Calculate floodplain characteristics
exp_bgs = phil_flood_record[phil_flood_record[!,string(2011)] .≥ 0.0,"GEOID"]

flpn_df = combine(groupby(price_df[price_df.GEOID .∈ Ref(exp_bgs),:],:time),
    Symbol.(shap_price_adata[6:8]) .=> (col -> sum(skipmissing(col))) .=> (string.(shap_price_adata[6:8])),
    Symbol.(shap_price_adata[9:11]) .=> (col -> sum(skipmissing(col))) .=> (string.(shap_price_adata[9:11])),
    Symbol.(shap_price_adata[15:17]) .=> (col -> sum(skipmissing(col))) .=> (string.(shap_price_adata[15:17]))
)
#Calculate average price
transform!(flpn_df,
    [:price_low, :cap_low] => ByRow(/) => :price_low_avg,
    [:price_med, :cap_med] => ByRow(/) => :price_med_avg,
    [:price_high, :cap_high] => ByRow(/) => :price_high_avg,
)
#normalize
for col in names(flpn_df)[2:end]
    flpn_df[!, col] = (flpn_df[!, col] .- flpn_df[1, col]) ./ flpn_df[1, col] .* 100
end

plot(flpn_df.time, [flpn_df.hh_low flpn_df.hh_high], label=["Low Income Pop." "High Income Pop."])
plot(flpn_df.time, [flpn_df.price_low_avg flpn_df.price_high_avg], label=["Low Inc. House Price" "High Inc. House Price"])




## Calculate relative demand over time
model = PhilPopABM(;input_params_dict...)
rel_demand = DataFrame()
t=1
while t < 21
    step!(model,1)
    bg_ids = [id for id in allids(model) if model[id] isa BlockGroup]
    n = length(bg_ids)

    bg_id = Vector{Int}(undef, n)  
    geo_id = Vector{Int}(undef, n)
    curr_util_low = Vector{Float64}(undef, n)
    curr_util_med = Vector{Float64}(undef, n)
    curr_util_high = Vector{Float64}(undef, n)
    dem_low = Vector{Float64}(undef, n)
    dem_med = Vector{Float64}(undef, n)
    dem_high = Vector{Float64}(undef, n)


    for (i, id) in enumerate(bg_ids)
        agent = model[id]
        bg_id[i] = id
        geo_id[i] = agent.GEOID
        curr_util_low[i] = agent.current_utility[1,1]
        curr_util_med[i] = agent.current_utility[2,2]
        curr_util_high[i] = agent.current_utility[3,3]
        dem_low[i] = agent.demand_exceeds_supply[1][model.tick]
        dem_med[i] = agent.demand_exceeds_supply[2][model.tick]
        dem_high[i] = agent.demand_exceeds_supply[3][model.tick]
    end

    dem_data = DataFrame(; bg_id, geo_id, curr_util_low,curr_util_med,curr_util_high,dem_low,dem_med,dem_high)
    dem_data.time = repeat([t],n)

    append!(rel_demand, dem_data)

    t += 1
    
end

exp_bgs = phil_flood_record[phil_flood_record[!,string(2011)] .> 0.0,"GEOID"]
rel_demand[rel_demand.geo_id .∈ Ref(exp_bgs),:]
#calculate avg utilities over time 
summary_df = combine(groupby(flpn_norm, :time),
    :hh_low => median => :hh_low_median,
    :hh_low => (x -> quantile(skipmissing(x), 0.025)) => :hh_low_lower,
    :hh_low => (x -> quantile(skipmissing(x), 0.975)) => :hh_low_upper,)

worst_bg_1981 = rand(phil_flood_record[phil_flood_record[!,string(1981)] .> 0.0,"GEOID"])[1]
phil_flood_record[phil_flood_record[!,"GEOID"] .== worst_bg_1981,"1981"]
rel_demand[rel_demand.geo_id .== worst_bg_1981,:]
price_df[price_df.GEOID .== worst_bg_1981,:][:,14:end]
model[geoid_to_bg[worst_bg_1981]].flood_hazard
model[geoid_to_bg[worst_bg_1981]].base_utility
model[geoid_to_bg[worst_bg_1981]].current_utility
##Evaluate price changes in select block group
worst_bg_2011 = rand(phil_flood_record[phil_flood_record[!,string(2011)] .> 0.2,"GEOID"])[1]

worst_bg_price = price_df[price_df.GEOID .== worst_bg_2011,:]
worst_bg_demand = rel_demand[rel_demand.geo_id .== worst_bg_2011, :]

##Look into demand calculation
model = PhilPopABM(;input_params_dict...)
step!(model,4)
#Update Year
model.tick += 1
#Reset outmigration count
map!(x->0, values(model.outmigrate))
#create new agents
CHANCE_C.AgentMigration(model; model.agent_creation...)
#Determine relocating HHAgents and potential moving locations
for id in collect(Agents.schedule(model))
    if model[id] isa HHAgent && model[id].bg_id < 1 #Dont involve HHAgents in Queues
        continue
    elseif model[id] isa CHANCE_C.Queue
        continue
    else
        CHANCE_C.agent_step!(model[id],model)
    end
end
CHANCE_C.LocationUpdate(model;grouped = model.build_develop[:grouped])
#run Housing Market to move HHAgents to desired locations

## Collect characteristics of all locations. Sort by utility value
loc_df = copy(model.df)
loc_df[!,:demand] = zeros(nrow(loc_df))
# Create a GEOID-to-BlockGroup lookup
geoid_to_bg = Dict{Int64, Int64}()
for bg in allagents(model)
    if bg isa BlockGroup
        geoid_to_bg[bg.GEOID] = bg.id
    end
end
##For each moving agent:
moving_agents = sort!([a for a in ids_in_position(model[0], model) if model[a] isa HHAgent], by=a -> model[a].income, rev=true)

sel_hh = []
sel_inc = []
sel_budg = []
sel_util = []
bg_util = []

stay_prob = 1.0
for ma in moving_agents
    #Subset to affordable and desirable locations
    bg_budget = subset(loc_df, :market_value => n -> n .<= model[ma].house_budget,
                :curr_utility => (c -> getindex.(c, model[ma].group) .>= first(values(model[ma].utility))),
                skipmissing=true, view = true
    )
    
    #Sort by utility value
    sort!(bg_budget, :curr_utility, by = x -> x[model[ma].group], rev = true)

    #check if selected block group is in bg_budget
    if length(bg_budget[bg_budget.GEOID .== worst_bg_2011,:GEOID]) > 0
        h_cat = bg_budget[bg_budget.GEOID .== worst_bg_2011,:income_cat][1]
        append!(sel_hh, model[ma].id)
        append!(sel_inc, model[ma].group)
        append!(sel_budg, model[ma].house_budget)
        append!(sel_util, first(values(model[ma].utility)))
        append!(bg_util, model[geoid_to_bg[worst_bg_2011]].current_utility[h_cat,model[ma].group])
    end

    #find first location with vacancy
    loc_ind = findfirst(x -> x > 0, bg_budget.available_units)
    #If there are no location options, agent moves back or outmigrates
    if isnothing(loc_ind)
        last_bg = model[first(keys(model[ma].utility))]
        if last_bg.id == -1
            remove_agent!(model[ma], model)
            continue
        elseif last_bg.available_units[model[ma].occ_cat] <= 0
            model.outmigrate[model[ma].group] += 1
            remove_agent!(model[ma], model)
            continue
        end
                    
        if rand(abmrng(model), Binomial(1, stay_prob)) == 1
            #Revert HHAgent Properties
            setproperty!(model[ma], :bg_id, last_bg.id)
            move_agent!(model[ma], last_bg.pos, model)
            #Update Last BG Properties
            last_bg.occupied_units[model[ma].occ_cat] += 1
            last_bg.available_units[model[ma].occ_cat] -= 1
            last_bg.population += getproperty(model[ma], :no_hhs_per_agent) * getproperty(model[ma], :hh_size)
            continue
        else
            model.outmigrate[model[ma].group] += 1
            remove_agent!(model[ma], model)
            continue
        end
    end
            
        
    #Get characteristics of block group
    last_bg_id = model[first(keys(model[ma].utility))].id
    new_bg_id = geoid_to_bg[bg_budget[loc_ind,:GEOID]]
    occ_cat = bg_budget[loc_ind,:income_cat]
    new_util = bg_budget[loc_ind, :curr_utility][model[ma].group] 
        
    #Move agent to new location
    move_agent!(model[ma], model[new_bg_id].pos, model)
    #Update bg_id, utility,  year of residence of agent
    setproperty!(model[ma], :bg_id, new_bg_id)
    setproperty!(model[ma], :occ_cat, occ_cat)
    setproperty!(model[ma], :utility, Dict(new_bg_id => new_util)) #Dict(bg_id => model[bg_id].current_utility[cat]))
    setproperty!(model[ma], :year_of_residence, model.tick)
    #update bg attributes
    model[new_bg_id].occupied_units[occ_cat] += 1
    model[new_bg_id].available_units[occ_cat] -= 1              
    model[new_bg_id].population += getproperty(model[ma],:no_hhs_per_agent) * getproperty(model[ma],:hh_size)
    #If moving agent is in-migrating, record in migrating agent dict
    if last_bg_id == -1
        model[new_bg_id].new_agents[model[ma].group] += 1
    end


    bg_budget[loc_ind, :available_units] -= 1
    #increase the interest count for selections at and above selected location
    bg_budget[1:loc_ind, :demand] .+= 1 
end

sel_df = DataFrame(;sel_hh, sel_inc, sel_budg, sel_util, bg_util)
model[geoid_to_bg[worst_bg_2011]].current_utility
#Update demand attributes for all BlockGroups 
for row in eachrow(loc_df)
    model[geoid_to_bg[row.GEOID]].demand_exceeds_supply[row.income_cat][model.tick] = row.demand - row.available_units
end

t_ma = rand(moving_agents)
model[t_ma].group
model[t_ma].house_budget
first(values(model[t_ma].utility))
bg_budget = subset(loc_df, :market_value => n -> n .<= model[t_ma].house_budget,
                :curr_utility => (c -> getindex.(c, model[t_ma].group) .>= first(values(model[t_ma].utility))),
                skipmissing=true, view = true
)
sort!(bg_budget, :curr_utility, by = x -> x[model[t_ma].group], rev = true)
bg_budget[:, [:GEOID,:market_value,:curr_utility, :available_units]]

loc_ind = findfirst(x -> x > 0, bg_budget.available_units)
bg_budget[1:loc_ind, :demand] .+= 1
#bg_budget_2[:, [:GEOID,:market_value,:curr_utility]]

getindex.(bg_budget.curr_utility, 2)