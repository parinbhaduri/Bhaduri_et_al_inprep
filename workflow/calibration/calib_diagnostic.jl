##File to ensure model runs sufficiently capture the output space across parameter combinations

#activate project environment
import Pkg
Pkg.activate(".")
Pkg.instantiate()

using Plots
using Plots.PlotMeasures
using CSV, DataFrames

using CairoMakie
using StatsPlots

#Read in Data
RUN_NO = 135842
par_combs = DataFrame(CSV.File(joinpath(@__DIR__,"data/param_comb_initial_$(RUN_NO).csv")))
calib_par_combs = DataFrame(CSV.File(joinpath(@__DIR__,"data/param_comb_final_$(RUN_NO)_mean_thresh_3.csv")))
sim_outputs = DataFrame(CSV.File(joinpath(@__DIR__,"data/calib_sim_output_$(RUN_NO).csv")))

phil_obs = DataFrame(CSV.File(joinpath(dirname(pwd()), "philadelphia-data","model_inputs", "calibration", "phil_obs_df.csv")))

param_cols = ["env_amen_l","price_inc_perc","rhea_coef","base_move",
"prop_l","prop_m","env_amen_m","risk_averse","build_inc_perc",
"env_amen_h","flood_coefficient","penalty","prop_h"]

### Look at top 10 quantile of highest avg variances
var_df_low = subset(par_combs, :var_l => x -> x .>= quantile(x, 0.9))
Dict(col => unique(var_df_low[!, col]) for col in param_cols)

var_df_med = subset(par_combs, :var_m => x -> x .>= quantile(x, 0.9))
Dict(col => unique(var_df_med[!, col]) for col in param_cols)

var_df_high = subset(par_combs, :var_h => x -> x .>= quantile(x, 0.9))
Dict(col => unique(var_df_high[!, col]) for col in param_cols)

##Look at the parameter combinations that have the highest low-inc and high-inc variance
comb_var = innerjoin(select(var_df_low, param_cols), select(var_df_high, param_cols), on = param_cols)
results = Dict()
for col in param_cols
    param_df = combine(groupby(comb_var, col), nrow => :count)
    results[col] = Dict(Pair.(param_df[!,col], param_df.count))
end
results

### Plot Distributions of Avg variances
stacked_df = stack(par_combs, [:var_l, :var_m, :var_h])
p = Plots.plot(size = (1000, 750), layout=(1,1), dpi = 300, plot_title = "Outputs by parameter")

StatsPlots.boxplot(p, repeat(["Med Income"], length(par_combs.var_m)), par_combs.var_m, outliers = true)
StatsPlots.boxplot(p, stacked_df.variable, stacked_df.value, outliers = false)




### Plot Distributions of outputs across different parameter combinations
function calib_spread(; output_df = sim_outputs, calib_df = calib_sim_outputs)
    p_low = Plots.plot(size = (1000, 750), layout=(2, 2), dpi = 300, plot_title = "Low Income Outputs")
    p_med = Plots.plot(size = (1000, 750), layout=(2, 2), dpi = 300, plot_title = "Middle Income Outputs")
    p_high = Plots.plot(size = (1000, 750), layout=(2, 2), dpi = 300, plot_title = "High Income Outputs")
    
    #param_val = unique(output_df[!, param_col])
    #df_low = subset(output_df, param_col => ByRow(isequal(param_val[1])))
    #df_high = subset(output_df, param_col => ByRow(isequal(param_val[2])))

    #For Low Income Populations
    histogram!(p_low[1], output_df.pop_prop_low, alpha = 0.5, label="Initial", fill = true,
    normalize = :density, legend_foreground_color = :transparent, left_margin = 5mm, bottom_margin = 5mm)
    histogram!(p_low[1], calib_df.pop_prop_low, alpha = 0.5, label = "Calibrated", fill = true, normalize = :density)
    vline!(p_low[1], phil_obs[!,"LOW_INC_PROP"], lw=3, label = "Observation")
    Plots.xlims!(p_low[1], 0.0, 0.5)
    Plots.ylabel!(p_low[1], "Count"; yguidefontsize=10)
    Plots.xlabel!(p_low[1], "Proportion of Low Inc. Population"; xguidefontsize=10)

    histogram!(p_low[2], output_df.pop_growth_low, alpha = 0.5, label="Initial", fill = true,
    normalize = :density, legend_foreground_color = :transparent, left_margin = 5mm, bottom_margin = 5mm)
    histogram!(p_low[2], calib_df.pop_growth_low, alpha = 0.5, label = "Calibrated", fill = true, normalize = :density)
    vline!(p_low[2], phil_obs[!,"LOW_INC_CHNG"], lw=3, label = "Observation")
    Plots.xlims!(p_low[2], -0.2, 0.3)
    Plots.ylabel!(p_low[2], "Count"; yguidefontsize=10)
    Plots.xlabel!(p_low[2], "Pop Growth of Low Inc. Population"; xguidefontsize=10)

    histogram!(p_low[3], output_df.moved_prop_low, alpha = 0.5, label="Initial", fill = true,
    normalize = :density, legend_foreground_color = :transparent, left_margin = 5mm, bottom_margin = 5mm)
    histogram!(p_low[3], calib_df.moved_prop_low, alpha = 0.5, label = "Calibrated", fill = true, normalize = :density)
    vline!(p_low[3], phil_obs[!,"LOW_SALE_PROP"], lw=3, label = "Observation")
    Plots.xlims!(p_low[3], 0, 0.6)
    Plots.ylabel!(p_low[3], "Count"; yguidefontsize=10)
    Plots.xlabel!(p_low[3], "Proportion of Low Inc. Housing Transactions"; xguidefontsize=10)

    histogram!(p_low[4], output_df.price_growth_low, alpha = 0.5, label="Initial", fill = true,
    normalize = :density, legend_foreground_color = :transparent, left_margin = 5mm, bottom_margin = 5mm)
    histogram!(p_low[4], calib_df.price_growth_low, alpha = 0.5, label = "Calibrated", fill = true, normalize = :density)
    vline!(p_low[4], phil_obs[!,"LOW_SALE_GROWTH"], lw=3, label = "Observation")
    Plots.xlims!(p_low[4], 0, 5.0)
    Plots.ylabel!(p_low[4], "Count"; yguidefontsize=10)
    Plots.xlabel!(p_low[4], "Price Growth of Low Inc. Housing"; xguidefontsize=10)

    #For Middle Income Populations
    histogram!(p_med[1], output_df.pop_prop_med, alpha = 0.5, label="Initial", fill = true,
    normalize = :density, legend_foreground_color = :transparent, left_margin = 5mm, bottom_margin = 5mm)
    histogram!(p_med[1], calib_df.pop_prop_med, alpha = 0.5, label = "Calibrated", fill = true, normalize = :density)
    vline!(p_med[1], phil_obs[!,"MED_INC_PROP"], lw=3, label = "Observation")
    Plots.xlims!(p_med[1], 0.0, 0.5)
    Plots.ylabel!(p_med[1], "Count"; yguidefontsize=10)
    Plots.xlabel!(p_med[1], "Proportion of Middle Inc. Population"; xguidefontsize=10)

    histogram!(p_med[2], output_df.pop_growth_med, alpha = 0.5, label="Initial", fill = true,
    normalize = :density, legend_foreground_color = :transparent, left_margin = 5mm, bottom_margin = 5mm)
    histogram!(p_med[2], calib_df.pop_growth_med, alpha = 0.5, label = "Calibrated", fill = true, normalize = :density)
    vline!(p_med[2], phil_obs[!,"MED_INC_CHNG"], lw=3, label = "Observation")
    Plots.xlims!(p_med[2], -0.2, 0.3)
    Plots.ylabel!(p_med[2], "Count"; yguidefontsize=10)
    Plots.xlabel!(p_med[2], "Pop Growth of Middle Inc. Population"; xguidefontsize=10)

    histogram!(p_med[3], output_df.moved_prop_med, alpha = 0.5, label="Initial", fill = true,
    normalize = :density, legend_foreground_color = :transparent, left_margin = 5mm, bottom_margin = 5mm)
    histogram!(p_med[3], calib_df.moved_prop_med, alpha = 0.5, label = "Calibrated", fill = true, normalize = :density)
    vline!(p_med[3], phil_obs[!,"MED_SALE_PROP"], lw=3, label = "Observation")
    Plots.xlims!(p_med[3], 0, 0.6)
    Plots.ylabel!(p_med[3], "Count"; yguidefontsize=10)
    Plots.xlabel!(p_med[3], "Proportion of Middle Inc. Housing Transactions"; xguidefontsize=10)

    histogram!(p_med[4], output_df.price_growth_med, alpha = 0.5, label="Initial", fill = true,
    normalize = :density, legend_foreground_color = :transparent, left_margin = 5mm, bottom_margin = 5mm)
    histogram!(p_med[4], calib_df.price_growth_med, alpha = 0.5, label = "Calibrated", fill = true, normalize = :density)
    vline!(p_med[4], phil_obs[!,"MED_SALE_GROWTH"], lw=3, label = "Observation")
    Plots.xlims!(p_med[4], 0, 5.0)
    Plots.ylabel!(p_med[4], "Count"; yguidefontsize=10)
    Plots.xlabel!(p_med[4], "Price Growth of Middle Inc. Housing"; xguidefontsize=10)

    #For High Income Populations
    histogram!(p_high[1], output_df.pop_prop_high, alpha = 0.5, label="Initial", fill = true,
    normalize = :density, legend_foreground_color = :transparent, left_margin = 5mm, bottom_margin = 5mm)
    histogram!(p_high[1], calib_df.pop_prop_high, alpha = 0.5, label = "Calibrated", fill = true, normalize = :density)
    vline!(p_high[1], phil_obs[!,"HIGH_INC_PROP"], lw=3, label = "Observation")
    Plots.xlims!(p_low[1], 0.0, 0.55)
    Plots.ylabel!(p_high[1], "Count"; yguidefontsize=10)
    Plots.xlabel!(p_high[1], "Proportion of High Inc. Population"; xguidefontsize=10)

    histogram!(p_high[2],output_df.pop_growth_high, alpha = 0.5, label="Initial", fill = true,
    normalize = :density, legend_foreground_color = :transparent, left_margin = 5mm, bottom_margin = 5mm)
    histogram!(p_high[2], calib_df.pop_growth_high, alpha = 0.5, label = "Calibrated", fill = true, normalize = :density)
    vline!(p_high[2], phil_obs[!,"HIGH_INC_CHNG"], lw=3, label = "Observation")
    Plots.xlims!(p_low[2], -0.6, 0.6)
    Plots.ylabel!(p_high[2], "Count"; yguidefontsize=10)
    Plots.xlabel!(p_high[2], "Pop Growth of High Inc. Population"; xguidefontsize=10)
    
    histogram!(p_high[3], output_df.moved_prop_high, alpha = 0.5, label="Initial", fill = true,
    normalize = :density, legend_foreground_color = :transparent, left_margin = 5mm, bottom_margin = 5mm)
    histogram!(p_high[3], calib_df.moved_prop_high, alpha = 0.5, label = "Calibrated", fill = true, normalize = :density)
    vline!(p_high[3], phil_obs[!,"HIGH_SALE_PROP"], lw=3, label = "Observation")
    Plots.xlims!(p_low[3], 0.0, 0.75)
    Plots.ylabel!(p_high[3], "Count"; yguidefontsize=10)
    Plots.xlabel!(p_high[3], "Proportion of High Inc. Housing Transactions"; xguidefontsize=10)

    histogram!(p_high[4],  output_df.price_growth_high, alpha = 0.5, label="Initial", fill = true,
    normalize = :density, legend_foreground_color = :transparent, left_margin = 5mm, bottom_margin = 5mm)
    histogram!(p_high[4],  calib_df.price_growth_high, alpha = 0.5, label = "Calibrated", fill = true, normalize = :density)
    vline!(p_high[4], phil_obs[!,"HIGH_SALE_GROWTH"], lw=3, label = "Observation")
    Plots.xlims!(p_low[4], 0, 5.0)
    Plots.ylabel!(p_high[4], "Count"; yguidefontsize=10)
    Plots.xlabel!(p_high[4], "Price Growth of High Inc. Housing"; xguidefontsize=10)
    
    return p_low, p_med, p_high
end


function output_spread(param_col; output_df = sim_outputs)
    p_low = Plots.plot(size = (1000, 750), layout=(2, 2), dpi = 300, plot_title = "Outputs by parameter : $(param_col)")
    p_med = Plots.plot(size = (1000, 750), layout=(2, 2), dpi = 300, plot_title = "Outputs by parameter : $(param_col)")
    p_high = Plots.plot(size = (1000, 750), layout=(2, 2), dpi = 300, plot_title = "Outputs by parameter : $(param_col)")
    
    param_val = unique(output_df[!, param_col])
    df_low = subset(output_df, param_col => ByRow(isequal(param_val[1])))
    df_high = subset(output_df, param_col => ByRow(isequal(param_val[2])))

    #For Low Income Populations
    histogram!(p_low[1], df_low.pop_prop_low, alpha = 0.5, label="Low Value", fill = true,
    normalize = :density, legend_foreground_color = :transparent, left_margin = 5mm, bottom_margin = 5mm)
    histogram!(p_low[1], df_high.pop_prop_low, alpha = 0.5, label = "High Value", fill = true, normalize = :density)
    vline!(p_low[1], phil_obs[!,"LOW_INC_PROP"], lw=3, label = "Observation")
    Plots.xlims!(p_low[1], 0.0, 0.5)
    Plots.ylabel!(p_low[1], "Count"; yguidefontsize=10)
    Plots.xlabel!(p_low[1], "Proportion of Low Inc. Population"; xguidefontsize=10)

    histogram!(p_low[2], df_low.pop_growth_low, alpha = 0.5, label="Low Value", fill = true,
    normalize = :density, legend_foreground_color = :transparent, left_margin = 5mm, bottom_margin = 5mm)
    histogram!(p_low[2], df_high.pop_growth_low, alpha = 0.5, label = "High Value", fill = true, normalize = :density)
    vline!(p_low[2], phil_obs[!,"LOW_INC_CHNG"], lw=3, label = "Observation")
    Plots.xlims!(p_low[2], -0.2, 0.3)
    Plots.ylabel!(p_low[2], "Count"; yguidefontsize=10)
    Plots.xlabel!(p_low[2], "Pop Growth of Low Inc. Population"; xguidefontsize=10)

    histogram!(p_low[3], df_low.moved_prop_low, alpha = 0.5, label="Low Value", fill = true,
    normalize = :density, legend_foreground_color = :transparent, left_margin = 5mm, bottom_margin = 5mm)
    histogram!(p_low[3], df_high.moved_prop_low, alpha = 0.5, label = "High Value", fill = true, normalize = :density)
    vline!(p_low[3], phil_obs[!,"LOW_SALE_PROP"], lw=3, label = "Observation")
    Plots.xlims!(p_low[3], 0, 0.6)
    Plots.ylabel!(p_low[3], "Count"; yguidefontsize=10)
    Plots.xlabel!(p_low[3], "Proportion of Low Inc. Housing Transactions"; xguidefontsize=10)

    histogram!(p_low[4], df_low.price_growth_low[df_low.price_growth_low .< 5], alpha = 0.5, label="Low Value", fill = true,
    normalize = :density, legend_foreground_color = :transparent, left_margin = 5mm, bottom_margin = 5mm)
    histogram!(p_low[4], df_high.price_growth_low[df_high.price_growth_low .< 5], alpha = 0.5, label = "High Value", fill = true, normalize = :density)
    vline!(p_low[4], phil_obs[!,"LOW_SALE_GROWTH"], lw=3, label = "Observation")
    Plots.xlims!(p_low[4], 0, 5.0)
    Plots.ylabel!(p_low[4], "Count"; yguidefontsize=10)
    Plots.xlabel!(p_low[4], "Price Growth of Low Inc. Housing"; xguidefontsize=10)

    #For Middle Income Populations
    histogram!(p_med[1], df_low.pop_prop_med, alpha = 0.5, label="Low Value", fill = true,
    normalize = :density, legend_foreground_color = :transparent, left_margin = 5mm, bottom_margin = 5mm)
    histogram!(p_med[1], df_high.pop_prop_med, alpha = 0.5, label = "High Value", fill = true, normalize = :density)
    vline!(p_med[1], phil_obs[!,"MED_INC_PROP"], lw=3, label = "Observation")
    Plots.xlims!(p_med[1], 0.0, 0.5)
    Plots.ylabel!(p_med[1], "Count"; yguidefontsize=10)
    Plots.xlabel!(p_med[1], "Proportion of Middle Inc. Population"; xguidefontsize=10)

    histogram!(p_med[2], df_low.pop_growth_med, alpha = 0.5, label="Low Value", fill = true,
    normalize = :density, legend_foreground_color = :transparent, left_margin = 5mm, bottom_margin = 5mm)
    histogram!(p_med[2], df_high.pop_growth_med, alpha = 0.5, label = "High Value", fill = true, normalize = :density)
    vline!(p_med[2], phil_obs[!,"MED_INC_CHNG"], lw=3, label = "Observation")
    Plots.xlims!(p_med[2], -0.2, 0.3)
    Plots.ylabel!(p_med[2], "Count"; yguidefontsize=10)
    Plots.xlabel!(p_med[2], "Pop Growth of Middle Inc. Population"; xguidefontsize=10)

    histogram!(p_med[3], df_low.moved_prop_med, alpha = 0.5, label="Low Value", fill = true,
    normalize = :density, legend_foreground_color = :transparent, left_margin = 5mm, bottom_margin = 5mm)
    histogram!(p_med[3], df_high.moved_prop_med, alpha = 0.5, label = "High Value", fill = true, normalize = :density)
    vline!(p_med[3], phil_obs[!,"MED_SALE_PROP"], lw=3, label = "Observation")
    Plots.xlims!(p_med[3], 0, 0.6)
    Plots.ylabel!(p_med[3], "Count"; yguidefontsize=10)
    Plots.xlabel!(p_med[3], "Proportion of Middle Inc. Housing Transactions"; xguidefontsize=10)

    histogram!(p_med[4], df_low.price_growth_med[df_low.price_growth_med .< 5], alpha = 0.5, label="Low Value", fill = true,
    normalize = :density, legend_foreground_color = :transparent, left_margin = 5mm, bottom_margin = 5mm)
    histogram!(p_med[4], df_high.price_growth_med[df_high.price_growth_med .< 5], alpha = 0.5, label = "High Value", fill = true, normalize = :density)
    vline!(p_med[4], phil_obs[!,"MED_SALE_GROWTH"], lw=3, label = "Observation")
    Plots.xlims!(p_med[4], 0, 5.0)
    Plots.ylabel!(p_med[4], "Count"; yguidefontsize=10)
    Plots.xlabel!(p_med[4], "Price Growth of Middle Inc. Housing"; xguidefontsize=10)

    #For High Income Populations
    histogram!(p_high[1], df_low.pop_prop_high, alpha = 0.5, label="Low Value", fill = true,
    normalize = :density, legend_foreground_color = :transparent, left_margin = 5mm, bottom_margin = 5mm)
    histogram!(p_high[1], df_high.pop_prop_high, alpha = 0.5, label = "High Value", fill = true, normalize = :density)
    vline!(p_high[1], phil_obs[!,"HIGH_INC_PROP"], lw=3, label = "Observation")
    Plots.xlims!(p_low[1], 0.0, 0.55)
    Plots.ylabel!(p_high[1], "Count"; yguidefontsize=10)
    Plots.xlabel!(p_high[1], "Proportion of High Inc. Population"; xguidefontsize=10)

    histogram!(p_high[2], df_low.pop_growth_high, alpha = 0.5, label="Low Value", fill = true,
    normalize = :density, legend_foreground_color = :transparent, left_margin = 5mm, bottom_margin = 5mm)
    histogram!(p_high[2], df_high.pop_growth_high, alpha = 0.5, label = "High Value", fill = true, normalize = :density)
    vline!(p_high[2], phil_obs[!,"HIGH_INC_CHNG"], lw=3, label = "Observation")
    Plots.xlims!(p_low[2], -0.6, 0.6)
    Plots.ylabel!(p_high[2], "Count"; yguidefontsize=10)
    Plots.xlabel!(p_high[2], "Pop Growth of High Inc. Population"; xguidefontsize=10)
    
    histogram!(p_high[3], df_low.moved_prop_high, alpha = 0.5, label="Low Value", fill = true,
    normalize = :density, legend_foreground_color = :transparent, left_margin = 5mm, bottom_margin = 5mm)
    histogram!(p_high[3], df_high.moved_prop_high, alpha = 0.5, label = "High Value", fill = true, normalize = :density)
    vline!(p_high[3], phil_obs[!,"HIGH_SALE_PROP"], lw=3, label = "Observation")
    Plots.xlims!(p_low[3], 0.0, 0.75)
    Plots.ylabel!(p_high[3], "Count"; yguidefontsize=10)
    Plots.xlabel!(p_high[3], "Proportion of High Inc. Housing Transactions"; xguidefontsize=10)

    histogram!(p_high[4],  df_low.price_growth_high[df_low.price_growth_high .< 5], alpha = 0.5, label="Low Value", fill = true,
    normalize = :density, legend_foreground_color = :transparent, left_margin = 5mm, bottom_margin = 5mm)
    histogram!(p_high[4],  df_high.price_growth_high[df_high.price_growth_high .< 5], alpha = 0.5, label = "High Value", fill = true, normalize = :density)
    vline!(p_high[4], phil_obs[!,"HIGH_SALE_GROWTH"], lw=3, label = "Observation")
    Plots.xlims!(p_low[4], 0, 5.0)
    Plots.ylabel!(p_high[4], "Count"; yguidefontsize=10)
    Plots.xlabel!(p_high[4], "Price Growth of High Inc. Housing"; xguidefontsize=10)
    
    return p_low, p_med, p_high
end


### Visualize Calibration Results ###
calib_sim_outputs = select(innerjoin(sim_outputs, calib_par_combs, on = param_cols), names(sim_outputs))
calib_low, calib_med, calib_high = calib_spread(;output_df = sim_outputs, calib_df = calib_sim_outputs)

for param_col in param_cols
    p_low, p_high = output_spread(param_col)
    savefig(p_low, joinpath(@__DIR__,"diagnostic_plots/output_spread_$(param_col)_low.png"))
    savefig(p_high, joinpath(@__DIR__,"diagnostic_plots/output_spread_$(param_col)_high.png"))
end

##Repeat for calibrated 
for param_col in param_cols
    p_low, p_med, p_high = output_spread(param_col;output_df=calib_sim_outputs)
    savefig(p_low, joinpath(@__DIR__,"diagnostic_plots/results_$(RUN_NO)/calib_output_spread_$(param_col)_low.png"))
    savefig(p_med, joinpath(@__DIR__,"diagnostic_plots/results_$(RUN_NO)/calib_output_spread_$(param_col)_high.png"))
    savefig(p_high, joinpath(@__DIR__,"diagnostic_plots/results_$(RUN_NO)/calib_output_spread_$(param_col)_high.png"))
end

param_val = unique(calib_sim_outputs[!, "prop_h"])