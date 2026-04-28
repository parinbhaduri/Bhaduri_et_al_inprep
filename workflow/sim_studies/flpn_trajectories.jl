import Pkg
Pkg.activate(".")
Pkg.instantiate()

using Plots, ColorSchemes
using ProgressMeter

### PARALLEL ENSEMBLE RUN ###
include(joinpath(dirname(@__DIR__), "src", "config_parallel.jl"))

param_path = joinpath(dirname(@__DIR__),"calibration","data/param_comb_final_mean_thresh_6_ens_250.csv")
calib_combs = DataFrame(CSV.File(param_path))[:,1:14]
shap_params = DataFrame(CSV.File(joinpath(dirname(@__DIR__),"shapley","data/shap_DESKTOP/param_runs_shap.csv")))#[:,1:14]
p_combs = collect((Tuple(row) for row in eachrow(calib_combs)))
output_params = collect(Symbol.(names(calib_combs)))

# Set up additional parameters
add_params = OrderedDict(
    :seed=>collect(range(1000,1249))
)

combs = [(p_combs[1]..., 2011, true,s) for s in values(add_params)...];
append!(output_params, [:flood_event_year, :flood_shock, :seed])
# Set up progress meter
progress = ProgressMeter.Progress(
    length(combs); 
    desc="Test Model Runs.. ",
    enabled=true,
    output=stderr,  # ProgressMeter outputs to stderr by default
    dt=5.0  # Update every 5 seconds
)

test_results = ProgressMeter.progress_pmap(combs; progress) do comb
    input_params_dict = Dict(output_params .=> comb)
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
    
end

#Combine model results into single df
adf = DataFrame()
seed = collect(1000:1249)
for (i,result) in enumerate(test_results)
    if nrow(result) > 0
        result.seed = repeat([seed[i]], nrow(result))
        append!(adf, result)   
    end
end
#Calculate floodplain characteristics
exp_bgs = phil_flood_record[phil_flood_record[!,string(2011)] .> 0.10,"GEOID"]

flpn_df = combine(groupby(adf[adf.GEOID .∈ Ref(exp_bgs),:],[:seed,:time]),
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
flpn_norm = transform(groupby(flpn_df,:seed),
    [Symbol(col) => (x -> (x .- x[1]) ./ x[1] .* 100) => Symbol(col) 
     for col in names(flpn_df)[3:end]]...
)

# First, calculate median and confidence intervals for each time point
summary_df = combine(groupby(flpn_norm, :time),
    :hh_low => median => :hh_low_median,
    :hh_low => (x -> quantile(skipmissing(x), 0.025)) => :hh_low_lower,
    :hh_low => (x -> quantile(skipmissing(x), 0.975)) => :hh_low_upper,
    :hh_med => median => :hh_med_median,
    :hh_med => (x -> quantile(skipmissing(x), 0.025)) => :hh_med_lower,
    :hh_med => (x -> quantile(skipmissing(x), 0.975)) => :hh_med_upper,
    :hh_high => median => :hh_high_median,
    :hh_high => (x -> quantile(skipmissing(x), 0.025)) => :hh_high_lower,
    :hh_high => (x -> quantile(skipmissing(x), 0.975)) => :hh_high_upper,

    :price_low_avg => median => :price_low_median,
    :price_low_avg => (x -> quantile(skipmissing(x), 0.025)) => :price_low_lower,
    :price_low_avg => (x -> quantile(skipmissing(x), 0.975)) => :price_low_upper,
    :price_med_avg => median => :price_med_median,
    :price_med_avg => (x -> quantile(skipmissing(x), 0.025)) => :price_med_lower,
    :price_med_avg => (x -> quantile(skipmissing(x), 0.975)) => :price_med_upper,
    :price_high_avg => median => :price_high_median,
    :price_high_avg => (x -> quantile(skipmissing(x), 0.025)) => :price_high_lower,
    :price_high_avg => (x -> quantile(skipmissing(x), 0.975)) => :price_high_upper,

)

p = Plots.plot(xlabel="Time", ylabel="Pop. Percent Change (%)", legend=:best)

# Column 1
Plots.plot!(summary_df.time, summary_df.hh_low_median, 
    ribbon=(summary_df.hh_low_median .- summary_df.hh_low_lower, 
            summary_df.hh_low_upper .- summary_df.hh_low_median),
    label="Low Income Pop.", linewidth=2, fillalpha=0.3, color=:blue)

# Column 2
Plots.plot!(summary_df.time, summary_df.hh_med_median, 
    ribbon=(summary_df.hh_med_median .- summary_df.hh_med_lower, 
            summary_df.hh_med_upper .- summary_df.hh_med_median),
    label="Middle Income Pop.", linewidth=2, fillalpha=0.3, color=:red)

# Column 3
Plots.plot!(summary_df.time, summary_df.hh_high_median, 
    ribbon=(summary_df.hh_high_median .- summary_df.hh_high_lower, 
            summary_df.hh_high_upper .- summary_df.hh_high_median),
    label="High Income Pop.", linewidth=2, fillalpha=0.3, color=:green)

display(p)



p2 = Plots.plot(xlabel="Time", ylabel="Price Percent Change (%)", legend=:best)

# Column 1
Plots.plot!(summary_df.time, summary_df.price_low_median, 
    ribbon=(summary_df.price_low_median .- summary_df.price_low_lower, 
            summary_df.price_low_upper .- summary_df.price_low_median),
    label="Low Val. Housing Price", linewidth=2, fillalpha=0.3, color=:blue)

# Column 2
Plots.plot!(summary_df.time, summary_df.price_med_median, 
    ribbon=(summary_df.price_med_median .- summary_df.price_med_lower, 
            summary_df.price_med_upper .- summary_df.price_med_median),
    label="Medium Val. Housing Price", linewidth=2, fillalpha=0.3, color=:red)

#Column 3
Plots.plot!(summary_df.time, summary_df.price_high_median, 
    ribbon=(summary_df.price_high_median .- summary_df.price_high_lower, 
            summary_df.price_high_upper .- summary_df.price_high_median),
    label="High Val. Housing Price", linewidth=2, fillalpha=0.3, color=:green)

display(p2)

