### Diagnostic Plots ###
# load environment and packages
import Pkg
Pkg.activate(".")
Pkg.instantiate()


using CSV, DataFrames
using Statistics 
using HDF5
using StatsBase


flpn_norm_low = DataFrame(CSV.File(joinpath(out_dir, "post_process","population","model_outcome_flpn_pop_norm_low.csv")))
flpn_norm_high = DataFrame(CSV.File(joinpath(out_dir, "post_process","population","model_outcome_flpn_pop_norm_high.csv")))
## Look at individual traces 
#Start with random runs
fig = Figure(size = (1000, 1000), fontsize = 24, pt_per_unit = 1, figure_padding = 18)
ax = Axis(fig[1, 1], ylabel = rich("Change in Population (%)";font=:bold), xlabel = rich("Time (years)";font=:bold), #xscale = log10,
titlesize = 28,  xgridvisible = false, titlealign = :center, title = "Change in Population in Exposed Areas")

random_row_indices = sample(axes(flpn_norm_low, 1), 10; replace = false)

for i in random_row_indices
    CairoMakie.lines!(ax, years, collect(flpn_norm_low[i,1:21]), color = Palette[4], alpha = 0.75, linewidth = 2)
    CairoMakie.lines!(ax, years, collect(flpn_norm_high[i,1:21]), color = Palette[2], alpha = 0.75, linewidth = 2)
end
#Create Legend
elem_1 = [LineElement(color = Palette[4], linestyle = :solid, linewidth = 5), PolyElement(color = (Palette[4], 0.35))]

elem_2 = [LineElement(color = Palette[2], linestyle = :solid, linewidth = 5), PolyElement(color = (Palette[2], 0.35))]

axislegend(ax, [elem_1, elem_2] , ["Low Income Population", "High Income Population"],
orientation = :vertical, framevisible = false, labelsize = 26, position = :lt)
display(fig)
#Grab runs with the greatest recovery 
#Grab runs with worst recovery 





##Look at trajectories for every population Distributions
model_params = DataFrame(CSV.File(joinpath(dirname(out_dir), "shap_DESKTOP","param_runs_shap.csv")))
pop_runs_dict = Dict(pop => findall(model_params.pop_no .== pop) for pop in collect(range(0,10)))
#Define output df
out_df = flpn_norm_high

fig = Figure(size = (1000, 1000), fontsize = 24, pt_per_unit = 1, figure_padding = 18)
ax = Axis(fig[1, 1], ylabel = rich("Change in Population (%)";font=:bold), xlabel = rich("Time (years)";font=:bold), #xscale = log10,
titlesize = 28,  xgridvisible = false, titlealign = :center, title = "Change in Population in Exposed Areas")

Palette = ColorSchemes.lipariS

for i in range(0,10,step=1)
    df_med = mapslices(x -> median(skipmissing(x)), Matrix(out_df[pop_runs_dict[i],:]), dims=1)
    df_quantiles = mapslices(x -> quantile(skipmissing(x), [0.025, 0.975]), Matrix(out_df[pop_runs_dict[i],:]), dims=1)

    #plot
    CairoMakie.lines!(ax, years, vec(df_med)[1:21], color = Palette[i+1], linewidth = 3.5)
    #, label = false)
    CairoMakie.band!(ax, years, df_quantiles[1,:][1:21], df_quantiles[2,:][1:21], color = (Palette[i+1], 0.35))
end

display(fig)

## Look at trajectories for worst and least hit block groups
include(joinpath(dirname(@__DIR__),"src","data_include.jl"))
##Read in data
out_dir = joinpath(@__DIR__,"data","shap_runs")
filename = "2011_abm_data_142772.h5"

h5file = h5open(joinpath(out_dir,filename), "r")
pop_dat = h5file["pop_data"]
#Grab the worst hit block groups
phil_flood_sorted = sort(phil_flood_record,"2011",rev=true)
ten_percent = ceil(Int, nrow(phil_flood_sorted) * 0.10) #equals 10% of all Block Groups
#Get the Block Group list for each category
worst_hit = phil_flood_sorted[1:ten_percent,"GEOID"]
least_hit = phil_flood_sorted[150:149+ten_percent,"GEOID"]
not_hit = phil_flood_record[phil_flood_record[!,"2011"] .== 0,"GEOID"]

#Select a exposure category
exp_hit = not_hit

#Translate to indices for array subsetting
exp_index = filter(x -> x !== nothing, indexin(exp_hit, h5file["GEOID"][:]))


exp_low = process_output(pop_dat;var_col=4, subset=true, index=exp_index)
init_exp_low = exp_low[:, 1]
exp_norm_low = (exp_low .- init_exp_low) ./ init_exp_low .* 100
exp_high = process_output(pop_dat;var_col=6, subset=true, index=exp_index)
init_exp_high = exp_high[:, 1]
exp_norm_high = (exp_high .- init_exp_high) ./ init_exp_high .* 100


low_med = mapslices(x -> median(skipmissing(x)), exp_norm_low, dims=1)
low_quantiles = mapslices(x -> quantile(skipmissing(x), [0.025, 0.975]), exp_norm_low, dims=1)

hi_med = mapslices(x -> median(skipmissing(x)), exp_norm_high, dims=1)
hi_quantiles = mapslices(x -> quantile(skipmissing(x), [0.025, 0.975]), exp_norm_high, dims=1)

fig = Figure(size = (1000,800), fontsize = 24, pt_per_unit = 1, figure_padding = 18)

ax1 = Axis(fig[1, 1], ylabel = rich("Change in Population (%)";font=:bold), xlabel = rich("Time (years)";font=:bold), #xscale = log10,
titlesize = 28,  xgridvisible = false, titlealign = :center, title = "Change in Population in Non- Exposed Areas") #xticks = ([10,100,1000], string.([10,100,1000])), limits = ((10,1000), nothing),

CairoMakie.lines!(ax1, years, vec(hi_med)[1:21], color = Palette[2], linewidth = 3.5)
#, label = false)
CairoMakie.band!(ax1, years, hi_quantiles[1,:][1:21], hi_quantiles[2,:][1:21], color = (Palette[2], 0.35))

CairoMakie.lines!(ax1, years, vec(low_med)[1:21], color = Palette[4], linewidth = 3.5)
#, label = false)
CairoMakie.band!(ax1, years, low_quantiles[1,:][1:21], low_quantiles[2,:][1:21], color = (Palette[4], 0.35))

#Create Legend
elem_1 = [LineElement(color = Palette[4], linestyle = :solid, linewidth = 5), PolyElement(color = (Palette[4], 0.35))]

elem_2 = [LineElement(color = Palette[2], linestyle = :solid, linewidth = 5), PolyElement(color = (Palette[2], 0.35))]

axislegend(ax1, [elem_1, elem_2] , ["Low Income Population", "High Income Population"],
orientation = :vertical, framevisible = false, labelsize = 26, position = :lt) #

display(fig)