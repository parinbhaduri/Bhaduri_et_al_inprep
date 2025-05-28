##File to ensure model runs sufficiently capture the output space across parameter combinations

#activate project environment
import Pkg
Pkg.activate(".")
Pkg.instantiate()

using Plots
using Plots.PlotMeasures
using CSV, DataFrames

#Read in Data
par_combs = DataFrame(CSV.File(joinpath(@__DIR__,"data/param_comb_initial.csv")))
sim_outputs = DataFrame(CSV.File(joinpath(@__DIR__,"data/calib_sim_output.csv")))

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

### Plot Distributions of outputs across different parameter combinations
function output_spread(param_col)
    p_low = Plots.plot(size = (1000, 750), layout=(2, 2), dpi = 300, plot_title = "Outputs by parameter : $(param_col)")
    p_high = Plots.plot(size = (1000, 750), layout=(2, 2), dpi = 300, plot_title = "Outputs by parameter : $(param_col)")
    
    param_val = unique(sim_outputs[!, param_col])
    df_low = subset(sim_outputs, param_col => ByRow(isequal(param_val[1])))
    df_high = subset(sim_outputs, param_col => ByRow(isequal(param_val[2])))

    #For Low Income Populations
    histogram!(p_low[1], df_low.pop_prop_low, alpha = 0.5, label="Low Value", fill = true,
    legend_foreground_color = :transparent, left_margin = 5mm, bottom_margin = 5mm)
    histogram!(p_low[1], df_high.pop_prop_low, alpha = 0.5, label = "High Value", fill = true)
    Plots.ylabel!(p_low[1], "Count"; yguidefontsize=10)
    Plots.xlabel!(p_low[1], "Proportion of Low Inc. Population"; xguidefontsize=10)

    histogram!(p_low[2], df_low.pop_growth_low, alpha = 0.5, label="Low Value", fill = true,
    legend_foreground_color = :transparent, left_margin = 5mm, bottom_margin = 5mm)
    histogram!(p_low[2], df_high.pop_growth_low, alpha = 0.5, label = "High Value", fill = true)
    Plots.ylabel!(p_low[2], "Count"; yguidefontsize=10)
    Plots.xlabel!(p_low[2], "Pop Growth of Low Inc. Population"; xguidefontsize=10)

    histogram!(p_low[3], df_low.moved_prop_low, alpha = 0.5, label="Low Value", fill = true,
    legend_foreground_color = :transparent, left_margin = 5mm, bottom_margin = 5mm)
    histogram!(p_low[3], df_high.moved_prop_low, alpha = 0.5, label = "High Value", fill = true)
    Plots.ylabel!(p_low[3], "Count"; yguidefontsize=10)
    Plots.xlabel!(p_low[3], "Proportion of Low Inc. Housing Transactions"; xguidefontsize=10)

    histogram!(p_low[4], df_low.price_growth_low[df_low.price_growth_low .< 5], alpha = 0.5, label="Low Value", fill = true,
    legend_foreground_color = :transparent, left_margin = 5mm, bottom_margin = 5mm)
    histogram!(p_low[4], df_high.price_growth_low[df_high.price_growth_low .< 5], alpha = 0.5, label = "High Value", fill = true)
    Plots.ylabel!(p_low[4], "Count"; yguidefontsize=10)
    Plots.xlabel!(p_low[4], "Price Growth of Low Inc. Housing"; xguidefontsize=10)

    #For High Income Populations
    histogram!(p_high[1], df_low.pop_prop_high, alpha = 0.5, label="Low Value", fill = true,
    legend_foreground_color = :transparent, left_margin = 5mm, bottom_margin = 5mm)
    histogram!(p_high[1], df_high.pop_prop_high, alpha = 0.5, label = "High Value", fill = true)
    Plots.ylabel!(p_high[1], "Count"; yguidefontsize=10)
    Plots.xlabel!(p_high[1], "Proportion of Low Inc. Population"; xguidefontsize=10)

    histogram!(p_high[2], df_low.pop_growth_high, alpha = 0.5, label="Low Value", fill = true,
    legend_foreground_color = :transparent, left_margin = 5mm, bottom_margin = 5mm)
    histogram!(p_high[2], df_high.pop_growth_high, alpha = 0.5, label = "High Value", fill = true)
    Plots.ylabel!(p_high[2], "Count"; yguidefontsize=10)
    Plots.xlabel!(p_high[2], "Pop Growth of Low Inc. Population"; xguidefontsize=10)

    histogram!(p_high[3], df_low.moved_prop_high, alpha = 0.5, label="Low Value", fill = true,
    legend_foreground_color = :transparent, left_margin = 5mm, bottom_margin = 5mm)
    histogram!(p_high[3], df_high.moved_prop_high, alpha = 0.5, label = "High Value", fill = true)
    Plots.ylabel!(p_high[3], "Count"; yguidefontsize=10)
    Plots.xlabel!(p_high[3], "Proportion of Low Inc. Housing Transactions"; xguidefontsize=10)

    histogram!(p_high[4],  df_low.price_growth_high[df_low.price_growth_high .< 5], alpha = 0.5, label="Low Value", fill = true,
    legend_foreground_color = :transparent, left_margin = 5mm, bottom_margin = 5mm)
    histogram!(p_high[4],  df_high.price_growth_high[df_high.price_growth_high .< 5], alpha = 0.5, label = "High Value", fill = true)
    Plots.ylabel!(p_high[4], "Count"; yguidefontsize=10)
    Plots.xlabel!(p_high[4], "Price Growth of Low Inc. Housing"; xguidefontsize=10)
    
    return p_low, p_high
end

for param_col in param_cols
    p_low, p_high = output_spread(param_col)
    savefig(p_low, joinpath(@__DIR__,"diagnostic_plots/output_spread_$(param_col)_low.png"))
    savefig(p_high, joinpath(@__DIR__,"diagnostic_plots/output_spread_$(param_col)_high.png"))
end

p_low, p_high = output_spread(param_cols[2])