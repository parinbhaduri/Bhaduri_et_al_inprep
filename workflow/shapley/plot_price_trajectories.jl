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

#Calculate median and 95% intervals
out_dir = joinpath(@__DIR__,"data","shap_runs")
fld_year = 2011
flpn_norm_low = DataFrame(CSV.File(joinpath(out_dir, "post_process","price","$(fld_year)_model_outcome_flpn_avg_price_norm_low.csv")))[:,1:21]
flpn_norm_high = DataFrame(CSV.File(joinpath(out_dir, "post_process","price","$(fld_year)_model_outcome_flpn_avg_price_norm_high.csv")))[:,1:21]

years = collect(range(0.0,20.0))
low_med = mapslices(x -> median(skipmissing(x)), Matrix(flpn_norm_low), dims=1)
low_quantiles = mapslices(x -> quantile(skipmissing(x), [0.025, 0.975]), Matrix(flpn_norm_low), dims=1)

hi_med = mapslices(x -> median(skipmissing(x)), Matrix(flpn_norm_high), dims=1)
hi_quantiles = mapslices(x -> quantile(skipmissing(x), [0.025, 0.975]), Matrix(flpn_norm_high), dims=1)

#Plot results
Palette = ColorSchemes.Hokusai3

fig = Figure(size = (1000,800), fontsize = 24, pt_per_unit = 1, figure_padding = 18)

ax1 = Axis(fig[1, 1], ylabel = rich("Change in Price (%)";font=:bold), xlabel = rich("Time (years)";font=:bold), #xscale = log10,
titlesize = 28,  xgridvisible = false, titlealign = :center, title = "Change in Housing Price in Exposed Areas") #xticks = ([10,100,1000], string.([10,100,1000])), limits = ((10,1000), nothing),

#Plot uncertainty intervals for population trajectories
CairoMakie.lines!(ax1, years, vec(hi_med)[1:21], color = Palette[2], linewidth = 3.5)
#, label = false)
CairoMakie.band!(ax1, years, hi_quantiles[1,:][1:21], hi_quantiles[2,:][1:21], color = (Palette[2], 0.25))

CairoMakie.lines!(ax1, years, vec(low_med)[1:21], color = Palette[4], linewidth = 3.5)
#, label = false)
CairoMakie.band!(ax1, years, low_quantiles[1,:][1:21], low_quantiles[2,:][1:21], color = (Palette[4], 0.25))

#Plot example traces
random_row_indices = sample(axes(flpn_norm_low, 1), 5; replace = false)

for i in random_row_indices
    CairoMakie.lines!(ax1, years, collect(flpn_norm_low[i,1:21]), color = Palette[4], linestyle =:dot, linewidth = 3)
    CairoMakie.lines!(ax1, years, collect(flpn_norm_high[i,1:21]), color = Palette[2], linestyle =:dot, linewidth = 3)
end
#CairoMakie.lines!(ax1, years, collect(flpn_norm_low[97012,1:21]), color = Palette[4], linestyle =:dot, linewidth = 4.5)
#CairoMakie.lines!(ax1, years, collect(flpn_norm_high[97012,1:21]), color = Palette[2], linestyle =:dot, linewidth = 4.5)

#Add Flood Shock
CairoMakie.vlines!(ax1, 4, color = :black, linecap = :round, linestyle = :dash, ymax = 0.8, linewidth = 2.5)
text!(ax1, 0, 440, text=rich("Flood Shock Occurrence", font = :italic), align = (:left, :center), fontsize = 22)

#Create Legend
elem_1 = [LineElement(color = Palette[4], linestyle = :solid, linewidth = 5), PolyElement(color = (Palette[4], 0.35))]
elem_2 = [LineElement(color = Palette[4], linestyle = :dot, linewidth = 5)]
elem_3 = [LineElement(color = Palette[2], linestyle = :solid, linewidth = 5), PolyElement(color = (Palette[2], 0.35))]
elem_4 = [LineElement(color = Palette[2], linestyle = :dot, linewidth = 5)]

axislegend(ax1, [elem_1, elem_2, elem_3, elem_4] , ["Low Income Price", "Low Income Trajectories", 
    "High Income Price", "High Income Trajectories"], nbanks = 2,
    orientation = :vertical, framevisible = false, labelsize = 26, position = :lt) #

display(fig)

CairoMakie.save(joinpath(pwd(),"figures", "shapley", "2011_avg_price_recovery.png"), fig)