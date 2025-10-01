# load environment and packages
import Pkg
Pkg.activate(".")
Pkg.instantiate()

using CSV, DataFrames # read CSV of Shapley indices
using DataFrames # data structure for indices
using Plots # plotting library
using ColorSchemes
using Measures # adjust margins with explicit measures
using Statistics # get mean function
using StatsPlots
#using CategoricalArrays


### Load Data ###
shap_indices = DataFrame(CSV.File(joinpath(out_dir, "post_process","shapley_indices","shap_indices_flpn_pop.csv")))
shap_ind_low = DataFrame(CSV.File(joinpath(out_dir, "post_process","shapley_indices","shap_indices_flpn_low_income.csv")))


function normalize_shap_groups(shap_ind)
    # assign parameters to groups by subsystem
    groups = Dict(
        "Parameters"    => ["env_amen_l","env_amen_m","env_amen_h",
                                "price_inc_perc","build_inc_perc","rhea_coef",
                                "base_move","risk_averse","penalty",
                                "prop_l","prop_m","prop_h",
                                "flood_mem","flood_coefficient"
                            ],
        "Population"    => ["pop_no"],
        "Internal Variability"   => ["seed"],
        #"Damage"   => ["damage_seed"],
        #"Flood Hazard"   => ["intensity", "freq"],
    
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
    group_order = ["Parameters","Population","Internal Variability"]
    shap_norm = shap_norm[indexin(group_order, shap_group.group), :]
    shap_permute = permutedims(shap_norm, 1)
    return shap_permute
end

shap_pop = normalize_shap_groups(shap_indices)
shap_pop_low = normalize_shap_groups(shap_ind_low)

# plot indices over time
group_colors = [:red,:blue,:green]

inch = 96
mmx = inch / 25.4
function plot_shapley(yrs,sdf,group_colors; leg=false)
    p_shap = areaplot(yrs, Matrix(sdf[!, Not(:group)]), xlabel="Year", ylabel="Relative Group Importance", 
            color_palette=group_colors, guidefontsize=16, tickfontsize=12, legend=leg, 
            label=permutedims(names(sdf)[2:end]), legendfontsize=12, fg_color_legend=false, 
            rightmargin=5mm, topmargin=5mm, leftmargin=5mm,
    )
    #annotate!(p_shap, 2030, 1.07, text("a", :left, 18, "Helvetica Bold"))
    Plots.xticks!(p_shap, 1980:5:2019)
    Plots.xlims!(p_shap, (1980, 2019))
    Plots.ylims!(p_shap, (0, 1))

    return p_shap
end

shap_pop_plot = plot_shapley(yrs,shap_pop,group_colors)
Plots.title!("Group Uncertainty Importance: Floodplain Pop.")
shap_low_plot = plot_shapley(yrs,shap_pop_low,group_colors;leg=:right)
Plots.title!("Group Uncertainty Importance: Low Income in Floodplain")


plt = Plots.plot(shap_pop_plot, shap_low_plot, layout=(2,1), dpi=300, size=(275mmx, 210mmx))

savefig(plt,joinpath(pwd(),"figures","shapley","shap_grid.png"))

##Plot for just One Year
shap_pop_2011 = shap_indices[:,"mean_2011"]
shap_share = shap_pop_2011 ./ sum(shap_pop_2011)

shap_2011 = Plots.bar(shap_indices.feature_name, shap_share,xrotation=45, label=false,
             xlabel="Features", ylabel="Relative Feature Importance",
)

Plots.title!("Feature Importance (Floodplain Pop.) in 2011")

savefig(shap_2011,joinpath(pwd(),"figures","shapley","shap_2011.png"))