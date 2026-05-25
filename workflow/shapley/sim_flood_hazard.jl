#activate project environment
import Pkg
Pkg.activate(".")
Pkg.instantiate()

using Distributed
addprocs(12, exeflags="--project=$(Base.active_project())")

@everywhere begin
    using CSV, DataFrames
    using Statistics
    using DataStructures
    using Agents
    using CHANCE_C
    using LinearAlgebra
end

### PARALLEL ENSEMBLE RUN ###
@everywhere begin
    include(joinpath(dirname(@__DIR__),"src/data_collect.jl"))
    include("shap_functions.jl")
end


### For changing flood hazard intensity
#Load flood hazard categories
haz_cat = DataFrame(CSV.File(joinpath(dirname(pwd()), "philadelphia-data/model_inputs", "phil_flood_hist_categories.csv")))

flood_years = vcat(2010,2013,1988) #1991,2018,1981,1996,1989,2011

fh_params = Dict(
    :seed=>1300,
    :flood_event_year=>flood_years,
    :risk_averse=>0.7, 
    :flood_mem=>10, 
    :build_inc_perc=>0.02,
    :price_inc_perc=>0.02, 
    :rhea_coef=>0.65,
    :base_move=>0.03, 
    :prop_l=>0.75,
    :env_amen_l=>0.25,
    :prop_m=>0.25,
    :env_amen_m=>0.75,
    :prop_h=>0.75,
    :env_amen_h=>0.25,
    :penalty=>0.75,
    :flood_coefficient=>0.25,
    
)

adf_intense,_ = paramscan(fh_params, PhilPopABM; parallel=true, showprogress=true, adata=shap_adata, n=19)

##Clean Data
queue_pos = adf_intense[adf_intense.agent_type .== Symbol("CHANCE_C.Queue"),:].pos[1:2]
subset!(adf_intense, :pos => x -> .!(x .∈ Ref(queue_pos)))
#Remove missing values
data_df =combine(groupby(adf_intense,[:flood_event_year, :time, :pos]),
    :id => minimum => :bg_id,
    :GEOID .=> (col -> minimum(skipmissing(col))) .=> :GEOID,
    Symbol.(shap_adata[3:end]) .=> (col -> sum(skipmissing(col))) .=> (string.(shap_adata[3:end]) .* "_sum"),
    #Symbol.(shap_adata[9:8]) .=> (col -> mean(skipmissing(col))) .=> (string.(shap_adata[6:end]) .* "_avg")
)

#=Subset to floodplain
phil_bg = DataFrame(CSV.File(joinpath(dirname(pwd()), "philadelphia-data/model_inputs", "phil_flood_bg_2019_nomiss_v1.csv")))
flpn_bgs = phil_bg[phil_bg.perc_flpn_area .> 0,:]
#Grab BGs impacted from the worst flood on record (2011)
bgs_2011 = phil_flood_record[phil_flood_record[!,"2011"] .> 0,:]
#Select block groups in floodplain that were exposed to the worst flood
flpn_bg_ids = unique(flpn_bgs[flpn_bgs.GEOID .∈ Ref(bgs_2011.GEOID),:].GEOID) #unique(flpn_bgs.GEOID)
=#
#flpn_index = indexin(flpn_bg_ids, h5file["GEOID"][:])
#Select floodplain-only BGs. Calculate statistics 
event_year = 2010
fdf = combine(groupby(data_df, ["flood_event_year","time"])) do group
    exp_bgs = phil_flood_record[phil_flood_record[!,string(event_year)] .> 0,"GEOID"]
    flpn = subset(group, "GEOID" => x -> (x .∈ Ref(exp_bgs)))
    (
        flpn_low = sum(flpn.hh_low_sum),
        flpn_med = sum(flpn.hh_med_sum),
        flpn_high = sum(flpn.hh_high_sum),
        flpn_price_low = sum(flpn.price_low_sum.*flpn.cap_low_sum)/sum(flpn.cap_low_sum),
        flpn_price_med = sum(flpn.price_med_sum.*flpn.cap_med_sum)/sum(flpn.cap_med_sum),
        flpn_price_high = sum(flpn.price_high_sum.*flpn.cap_high_sum)/sum(flpn.cap_high_sum),
    )
end

##Plot results
#create subsets for hazard categories
#major = flood_years[[3,6,9]]
#medium = flood_years[[2,5,8]]
minor = flood_years[[1,2,3]]

function flpn_pop_plot(df, year_set; leg = :outertopright, color = palette(:berlin10), lim = (0,12000), p_lim = (1e5,1e6))
    ddf = subset(df, "flood_event_year" => x -> (x .∈ Ref(year_set)))
    ddf.flood_event_year = replace(ddf.flood_event_year, year_set[1] => "min", year_set[2] => "med", year_set[3] => "max") 

    plots = []
    push!(plots, Plots.plot(ddf.time, ddf.flpn_low, group = ddf[:, "flood_event_year"], legend=leg, palette = color, linewidth=2, ylimits=lim, title="Low Income Pop."))
    push!(plots, Plots.plot(ddf.time, ddf.flpn_price_low, group = ddf[:, "flood_event_year"], legend=leg, palette = color, linewidth=2, ylimits=p_lim, title="Avg. Housing Price (low)"))

    push!(plots, Plots.plot(ddf.time, ddf.flpn_med, group = ddf[:, "flood_event_year"], legend=leg, palette = color, linewidth=2, ylimits=lim, title="Medium Income Pop."))
    push!(plots, Plots.plot(ddf.time, ddf.flpn_price_med, group = ddf[:, "flood_event_year"], legend=leg, palette = color, linewidth=2, ylimits=p_lim, title="Avg. Housing Price (med)"))

    push!(plots, Plots.plot(ddf.time, ddf.flpn_high, group = ddf[:, "flood_event_year"], legend=leg, palette = color, linewidth=2, ylimits=lim, title="High Income Pop."))
    push!(plots, Plots.plot(ddf.time, ddf.flpn_price_high, group = ddf[:, "flood_event_year"], legend=leg, palette = color, linewidth=2, ylimits=p_lim, title="Avg. Housing Price (high)"))
    
    return plots
end



using Plots
using ColorSchemes

dat_df = copy(fdf)

maj_pop_plots = flpn_pop_plot(dat_df,major)
plot(maj_pop_plots..., layout=(3, 2), size = (1100, 1000), plot_title = "Flpn. Dynamics with Major Flood Shock")

med_pop_plots = flpn_pop_plot(dat_df,medium)
plot(med_pop_plots..., layout=(3, 2), size = (1100, 1000), plot_title = "Flpn. Dynamics with Medium Flood Shock")

min_pop_plots = flpn_pop_plot(dat_df,minor)
Plots.plot(min_pop_plots..., layout=(3, 2), size = (1100, 1000), plot_title = "Flpn. Dynamics with Minor Flood Shock")