# load environment and packages
import Pkg
Pkg.activate(".")
Pkg.instantiate()


using CSV, DataFrames
using Statistics 
using HDF5
using ProgressMeter
using StatsBase
using CairoMakie
using ColorSchemes
using Random

Random.seed!(1)
#Subset data by param combination
param_path = joinpath(dirname(@__DIR__),"calibration","data/param_comb_final_mean_thresh_6_ens_250.csv")
calib_combs = DataFrame(CSV.File(param_path))[:,1:14]
calib_combs.pop_no .= [0]
shap_params = DataFrame(CSV.File(joinpath(dirname(@__DIR__),"shapley","data/shap_DESKTOP/param_runs_shap.csv")))


#Calculate median and 95% intervals
out_dir = joinpath(@__DIR__,"data","shap_runs")
fld_year = 2013
fld_event = "Low"
flpn_norm_low = DataFrame(CSV.File(joinpath(out_dir, "post_process","population",fld_event,"$(fld_year)_model_outcome_flpn_pop_norm_low.csv")))[:,1:21]
#flpn_norm_med = DataFrame(CSV.File(joinpath(out_dir, "post_process","population","High","$(fld_year)_model_outcome_flpn_pop_norm_med.csv")))[:,1:21]
flpn_norm_high = DataFrame(CSV.File(joinpath(out_dir, "post_process","population",fld_event,"$(fld_year)_model_outcome_flpn_pop_norm_high.csv")))[:,1:21]
flpn_norm_tot = DataFrame(CSV.File(joinpath(out_dir, "post_process","population",fld_event,"$(fld_year)_model_outcome_flpn_pop_norm.csv")))[:,1:21]

#Load counterfactual data
counter_low = DataFrame(CSV.File(joinpath(out_dir, "post_process","population","High","1982_model_outcome_flpn_pop_norm_low.csv")))
counter_high = DataFrame(CSV.File(joinpath(out_dir, "post_process","population","High","1982_model_outcome_flpn_pop_norm_high.csv")))
counter_tot = DataFrame(CSV.File(joinpath(out_dir, "post_process","population","High","1982_model_outcome_flpn_pop_norm.csv")))

#Count no. of realizations where pop recovers and does not recover for low and high Income
(sum(flpn_norm_low[:,end] .>= 0)/length(flpn_norm_low[:,end])) * 100
(sum(flpn_norm_high[:,end] .>= 0)/length(flpn_norm_high[:,end])) * 100

#Calculate avg % diff in pop bw shock scenario and counterfactual
println("Avg. % difference in population between $(fld_event)($(fld_year)) flood event and counterfactual scenario:", mean(flpn_norm_tot[:,end]) - mean(counter_tot[:,end]))

#High income pop
println("High Income Population quick stats")
println("Average pop. growth: $(fld_year) - ", mean(flpn_norm_high[:,end]), "%")
println("Average pop. growth: Counterfactual - ", mean(counter_high[:,end]), "%")
#Total pop
println("Total Population quick stats")
println("Average pop. growth: $(fld_year) - ", mean(flpn_norm_tot[:,end]), "%")
println("Average pop. growth: Counterfactual - ", mean(counter_tot[:,end]), "%")
(mean(flpn_norm_tot[:,end]) - mean(counter_tot[:,end]))


years = collect(range(0.0,20.0))
low_med = mapslices(x -> median(skipmissing(x)), Matrix(flpn_norm_low), dims=1)
low_quantiles = mapslices(x -> quantile(skipmissing(x), [0.025, 0.975]), Matrix(flpn_norm_low), dims=1)

#med_med = mapslices(x -> median(skipmissing(x)), Matrix(flpn_norm_med), dims=1)
#med_quantiles = mapslices(x -> quantile(skipmissing(x), [0.025, 0.975]), Matrix(flpn_norm_med), dims=1)

hi_med = mapslices(x -> median(skipmissing(x)), Matrix(flpn_norm_high), dims=1)
hi_quantiles = mapslices(x -> quantile(skipmissing(x), [0.025, 0.975]), Matrix(flpn_norm_high), dims=1)

#Plot results
Palette = ColorSchemes.Hokusai3
hc = Palette[2]#High Income Color
lc = Palette[4]#low Income Color

fig = Figure(size = (1000,750), fontsize = 24, pt_per_unit = 1, figure_padding = 18)
ga = fig[1, 1] = GridLayout()
gb = fig[1, 2] = GridLayout()
ax1 = Axis(ga[1, 1], ylabel = rich("Change in Population (%)";font=:bold), xlabel = rich("Time (years)";font=:bold), #xscale = log10,
titlesize = 28,  xgridvisible = false) #, titlealign = :center, title = "Change in Population in Exposed Areas" #xticks = ([10,100,1000], string.([10,100,1000])), limits = ((10,1000), nothing),
CairoMakie.xlims!(ax1, low = 0, high=20)

ax2 = Axis(gb[1, 1], xticks = [1,2], xgridvisible = false)  #xlabel = rich("Population Group"; font = :bold),
hideydecorations!(ax2, grid = false)
linkyaxes!(ax1, ax2)

#Plot uncertainty intervals for population trajectories
CairoMakie.lines!(ax1, years, vec(hi_med)[1:21], color = hc, linewidth = 3.5)
#, label = false)
CairoMakie.band!(ax1, years, hi_quantiles[1,:][1:21], hi_quantiles[2,:][1:21], color = ( hc, 0.15))

#CairoMakie.lines!(ax1, years, vec(med_med)[1:21], color = Palette[3], linewidth = 3.5)
#, label = false)
#CairoMakie.band!(ax1, years, med_quantiles[1,:][1:21], med_quantiles[2,:][1:21], color = (Palette[3], 0.25))

CairoMakie.lines!(ax1, years, vec(low_med)[1:21], color = lc, linewidth = 3.5)
#, label = false)
CairoMakie.band!(ax1, years, low_quantiles[1,:][1:21], low_quantiles[2,:][1:21], color = (lc, 0.15))

#Plot example traces

random_row_indices = sample(axes(flpn_norm_low, 1), 5; replace = false)

for i in random_row_indices
    CairoMakie.lines!(ax1, years, collect(flpn_norm_low[i,1:21]), color = lc, linestyle =:dot, linewidth = 2.0)
    CairoMakie.lines!(ax1, years, collect(flpn_norm_high[i,1:21]), color =  hc, linestyle =:dot, linewidth = 2.0)
end
#CairoMakie.lines!(ax1, years, collect(flpn_norm_low[97012,1:21]), color = Palette[4], linestyle =:dot, linewidth = 4.5)
#CairoMakie.lines!(ax1, years, collect(flpn_norm_high[97012,1:21]), color = Palette[2], linestyle =:dot, linewidth = 4.5)
#Add Horiz Line at 0
CairoMakie.hlines!(ax1, 0, color = :black, linecap = :round, xmax = 20, linewidth = 2.5)
#Add Flood Shock
CairoMakie.vlines!(ax1, 4, color = :black, linecap = :round, linestyle = :dash, ymax = 0.8, linewidth = 2.5)
text!(ax1, 4.5, 52, text=rich("Flood Shock Occurrence", font = :italic), align = (:left, :center), fontsize = 22)

#Create Legend
elem_1 = [PolyElement(color = (lc, 0.5))]
elem_2 = [PolyElement(color = (hc, 0.5))]
elem_3 = [LineElement(color = :black, linestyle = :solid, linewidth = 5)]
elem_4 = [LineElement(color =  :black, linestyle = :dot, linewidth = 5)]

axislegend(ax1, [elem_1, elem_2,  elem_3, elem_4] , ["Low Income", "High Income",  
        "Median Trajectory", "Individual Trajectory"], nbanks = 2,
    orientation = :vertical, framevisible = false, labelsize = 26, position = :lt
) #

#plot Boxplot of final pop_change range (w/ counterfactual)
flpn_norm_low_final = flpn_norm_low[:,end]
flpn_norm_high_final = flpn_norm_high[:,end]
counter_low_final = counter_low[:,end]
counter_high_final = counter_high[:,end]
data_len = length(flpn_norm_low_final)


data = vcat(flpn_norm_low_final, flpn_norm_high_final, counter_low_final, counter_high_final)
groups = repeat([1,2],inner=data_len*2)
dodges = repeat([1,2,1,2],inner=data_len)
colors = vcat(
    fill(lc, data_len),
    fill(hc, data_len),
    fill((lc, 0.5), data_len),
    fill((hc, 0.35), data_len)
);

CairoMakie.boxplot!(ax2, groups, data, 
    dodge = dodges,
    color = colors,
    strokewidth = 1,
    width = 0.5,  # Make boxes narrower
    gap = 0.05,
)

ax2.xticks = (1:2, ["Final Year", "Final Year\n(No Flood)"])

colsize!(fig.layout, 1, Auto(3))
colgap!(fig.layout, 1, 40)
colsize!(fig.layout, 2, Auto(1))

display(fig)


CairoMakie.save(joinpath(pwd(),"figures", "recovery", "$(fld_year)_pop_recovery_w_box.png"), fig)
