# load environment and packages
import Pkg
Pkg.activate(".")
Pkg.instantiate()

using CSV, DataFrames # read CSV of Shapley indices
using Plots # plotting library
using ColorSchemes
using Measures # adjust margins with explicit measures
using Statistics # get mean function
using StatsPlots
#using CategoricalArrays


### Load Data ###
out_dir = joinpath(@__DIR__,"data","shap_runs")
shap_ind_high = DataFrame(CSV.File(joinpath(out_dir, "post_process","shapley_indices","test_2011_shap_indices_flpn_pop_norm_high.csv")))
shap_ind_low = DataFrame(CSV.File(joinpath(out_dir, "post_process","shapley_indices","test_2011_shap_indices_flpn_pop_norm_low.csv")))


function normalize_shap_groups(shap_ind)
    # assign parameters to groups by subsystem
    groups = Dict(
        "Flood Risk Perception"    => ["risk_averse", "flood_mem","flood_coefficient"],
        "Housing Market"    => ["price_inc_perc","build_inc_perc","rhea_coef","base_move"],                    
        "Location Preference"    => ["env_amen_l","env_amen_m","env_amen_h",
                                "penalty", "prop_l","prop_m","prop_h",
                            ],
        "Population"    => ["pop_no"],
        "Internal Stochasticity"   => ["seed"],
        #"Damage"   => ["damage_seed"],
        "Flood Hazard"   => ["fld_extents"],
    
    )

    group_col = Array{String}(undef, size(shap_ind)[1]) # preallocate space
    for (key,value) in groups # loop through dictionary and create vector with group labels
        indices = findall((in)(value), shap_ind.feature_name) # find indices for current group
        group_col[indices] .= key # fill the indices with the name of current group
    end
    shap_ind.group = group_col # add group column to dataframe   

    # sum shapley effects by group for each year
    shap_group = groupby(shap_ind[!, Not(:feature_name)], :group)
    shap_group = combine(shap_group, Not(:group) .=> sum)
    # normalize so the grouped shapley sums equal 1
    shap_norm = mapcols(x -> x / sum(x), shap_group[!, Not([:group])])
    insertcols!(shap_norm, 1, :group => shap_group.group)
    group_order = ["Flood Hazard","Flood Risk Perception","Housing Market", "Location Preference","Population","Internal Stochasticity"]
    shap_norm = shap_norm[indexin(group_order, shap_group.group), :]
    shap_permute = permutedims(shap_norm, 1)
    return shap_permute
end

shap_pop = normalize_shap_groups(shap_indices)
shap_pop_high = normalize_shap_groups(shap_ind_high)
shap_pop_low = normalize_shap_groups(shap_ind_low)

# plot indices over time
group_colors = ColorSchemes.mk_8[2:end]

inch = 96
mmx = inch / 25.4
function plot_shapley(yrs,sdf,group_colors; leg=false)
    p_shap = areaplot(yrs, Matrix(sdf[!, Not(:group)]), xlabel="Year", ylabel="Relative Group Importance", 
            color_palette=group_colors, guidefontsize=16, tickfontsize=12, legend=leg, legendcolumns=2,
            label=permutedims(names(sdf)[2:end]), legendfontsize=12, fg_color_legend=false, 
            rightmargin=5mm, topmargin=5mm, leftmargin=5mm,
    )
    #annotate!(p_shap, 2030, 1.07, text("a", :left, 18, "Helvetica Bold"))
    Plots.vline!(p_shap, [4], ls=:dash, lc=:black, lw=3, label = "Flood Shock")
    Plots.xticks!(p_shap, yrs[1]:5:yrs[end])
    Plots.xlims!(p_shap, (yrs[1],yrs[end]))
    Plots.ylims!(p_shap, (0, 1))

    return p_shap
end

yrs = collect(range(0,20))
shap_hi_plot = plot_shapley(yrs,shap_pop_high,group_colors; leg=:outerbottom)
#Plots.title!("Group Uncertainty Importance: Floodplain Pop.")
savefig(shap_hi_plot,joinpath(pwd(),"figures","shapley","test_2011_shap_grid_flpn_pop_norm_high_w_leg.png"))

shap_low_plot = plot_shapley(yrs,shap_pop_low,group_colors; leg=false)
#Plots.title!("Group Uncertainty Importance: Low Income in Floodplain")
savefig(shap_low_plot,joinpath(pwd(),"figures","shapley","test_2011_shap_grid_flpn_pop_norm_low.png"))

#plt = Plots.plot(shap_hi_plot, shap_low_plot,  layout = (1,2), dpi=300, size=(275mmx, 210mmx))

savefig(plt,joinpath(pwd(),"figures","shapley","shap_grid.png"))

##Plot for just One Year
shap_pop_2011 = shap_indices[:,"mean_2011"]
shap_share = shap_pop_2011 ./ sum(shap_pop_2011)

shap_2011 = Plots.bar(shap_indices.feature_name, shap_share,xrotation=45, label=false,
             xlabel="Features", ylabel="Relative Feature Importance",
)

Plots.title!("Feature Importance (Floodplain Pop.) in 2011")

savefig(shap_2011,joinpath(pwd(),"figures","shapley","shap_2011.png"))