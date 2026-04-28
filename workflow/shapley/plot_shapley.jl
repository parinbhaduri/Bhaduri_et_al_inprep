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
out_dir = joinpath(@__DIR__,"data","shap_runs") #"shapley_indices",

# assign parameters to groups by subsystem
out_groups = Dict(
        "Flood Memory and Aversion"    => ["risk_averse", "flood_mem","flood_coefficient"],
        "Housing Market"    => ["price_inc_perc","build_inc_perc","rhea_coef","base_move"],                    
        "Location Amenities"    => ["env_amen_l","env_amen_m","env_amen_h",
                                "penalty", "prop_l","prop_m","prop_h",
                            ],
        "Initial Population Distribution"    => ["pop_no"],
        "Internal Stochasticity"   => ["seed"],
        "Flood Hazard"   => ["fld_extents"],    
)
out_order = ["Flood Hazard","Flood Memory and Aversion", "Location Amenities", "Housing Market", "Initial Population Distribution","Internal Stochasticity"]

function normalize_shap_groups(shap_ind; normalize=true, groups = out_groups, order = out_order)
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
    if normalize
        shap_norm = mapcols(x -> x / sum(x), shap_group[!, Not([:group])])
    else
        shap_norm = mapcols(x -> x, shap_group[!, Not([:group])])
    end

    insertcols!(shap_norm, 1, :group => shap_group.group)
    
    shap_norm = shap_norm[indexin(order, shap_group.group), :]
    shap_permute = permutedims(shap_norm, 1)
    return shap_permute
end

function plot_shapley(yrs,sdf,group_colors; leg=false, normalize=true)
    if normalize
        y_label = "Relative Group Importance"
        y_lims = (0,1)
    else
        y_label = "Absolute Group Importance"
        y_lims = (0,100)
    end
    p_shap = areaplot(yrs, Matrix(sdf[!, Not(:group)]), xlabel="Year", ylabel=y_label, 
            color_palette=group_colors, guidefontsize=16, tickfontsize=12, legend=leg, legendcolumns=2,
            label=permutedims(names(sdf)[2:end]), legendfontsize=12, fg_color_legend=false, 
            rightmargin=5mm, topmargin=5mm, leftmargin=5mm,
    )
    #annotate!(p_shap, 2030, 1.07, text("a", :left, 18, "Helvetica Bold"))
    Plots.vline!(p_shap, [4], ls=:dash, lc=:black, lw=3, label = "Flood Shock")
    Plots.xticks!(p_shap, yrs[1]:5:yrs[end])
    Plots.xlims!(p_shap, (yrs[1],yrs[end]))
    Plots.ylims!(p_shap, y_lims)

    return p_shap
end

#shap_pop = normalize_shap_groups(shap_indices)
# plot indices over time
group_colors = vcat(ColorSchemes.mk_8[2:3],ColorSchemes.mk_8[5],ColorSchemes.mk_8[4],ColorSchemes.mk_8[6:end])

inch = 96
mmx = inch / 25.4

yrs = collect(range(0,20))
#Get stats

for (outcome,out_name) in zip(["population","price"],["pop","price"])
    for fld_event in ["High", "Medium","Low"] #
        shap_ind_high = DataFrame(CSV.File(joinpath(out_dir, "post_process","shapley_indices","$(fld_event)_fld_shap_indices_flpn_$(out_name)_norm_high.csv")))
        shap_ind_low = DataFrame(CSV.File(joinpath(out_dir, "post_process","shapley_indices","$(fld_event)_fld_shap_indices_flpn_$(out_name)_norm_low.csv")))
        shap_out_high = normalize_shap_groups(shap_ind_high)
        shap_out_low = normalize_shap_groups(shap_ind_low)
        shap_out_high_abs = normalize_shap_groups(shap_ind_high; normalize=false)
        shap_out_low_abs = normalize_shap_groups(shap_ind_low; normalize=false)
        if outcome == "population"
            println("Relative Flood Mem. Contribution for High Inc. Trends ($(fld_event)):",maximum(filter(!isnan,shap_out_high[:,"Flood Memory and Aversion"])))
            println("Relative Flood Mem. Contribution for Low Inc. Trends ($(fld_event)):",maximum(filter(!isnan,shap_out_low[:,"Flood Memory and Aversion"])))
        end

        #=
        shap_hi_plot = plot_shapley(yrs,shap_out_high,group_colors; leg=false)
        #Plots.title!("Group Uncertainty Importance: Floodplain Pop.")
        savefig(shap_hi_plot,joinpath(pwd(),"figures","shapley",outcome,"$(fld_event)_fld_shap_grid_flpn_$(out_name)_norm_high_rel.png"))

        shap_low_plot = plot_shapley(yrs,shap_out_low,group_colors; leg=false)
        #Plots.title!("Group Uncertainty Importance: Low Income in Floodplain")
        savefig(shap_low_plot,joinpath(pwd(),"figures","shapley",outcome,"$(fld_event)_fld_shap_grid_flpn_$(out_name)_norm_low_rel.png"))
        
        ## Plot Absolute Importance 
        shap_hi_plot = plot_shapley(yrs,shap_out_high_abs,group_colors; leg=false, normalize=false)
        #Plots.title!("Group Uncertainty Importance: Floodplain Pop.")
        savefig(shap_hi_plot,joinpath(pwd(),"figures","shapley",outcome,"$(fld_event)_fld_shap_grid_flpn_$(out_name)_norm_high_abs.png"))

        shap_low_plot = plot_shapley(yrs,shap_out_low_abs,group_colors; leg=false, normalize=false)
        #Plots.title!("Group Uncertainty Importance: Low Income in Floodplain")
        savefig(shap_low_plot,joinpath(pwd(),"figures","shapley",outcome,"$(fld_event)_fld_shap_grid_flpn_$(out_name)_norm_low_abs.png"))
        =#
    end
end
img_files = filter(file -> occursin(Regex("rel.png"),file), readdir(joinpath(pwd(),"figures","shapley","population")))
[FileIO.load(joinpath(pwd(),"figures","shapley","population",img)) for img in img_files]
### Create Multi-Panel Plot of Shapley evolution 
using CairoMakie
using FileIO
for outcome in ["population","price"]
    for scale in ["rel","abs"]

        img_files = filter(file -> occursin(Regex("$(scale).png"),file), readdir(joinpath(pwd(),"figures","shapley",outcome)))
        # Move Minor Events to end of file list
        min_imgs = [img_files[3], img_files[4]]
        deleteat!(img_files,(3,4))
        append!(img_files, min_imgs)
        # create Figure 
        fig = Figure(size = (1000, 900), fontsize = 20, pt_per_unit = 1, figure_padding = 18)

        # Add column labels at the top (row 1)
        if outcome == "population"
            Label(fig[1, 2], "Low Income Population", fontsize = 22, font = :bold)
            Label(fig[1, 3], "High Income Population", fontsize = 22, font = :bold)
        else
            Label(fig[1, 2], "Low Value Price", fontsize = 22, font = :bold)
            Label(fig[1, 3], "High Value Price", fontsize = 22, font = :bold)
        end
        

        # Add row labels on the left (column 1) 
        Label(fig[2, 1], "Extreme", fontsize = 22, font = :bold, rotation = π/2)
        Label(fig[3, 1], "Moderate", fontsize = 22, font = :bold, rotation = π/2)
        Label(fig[4, 1], "Nuisance", fontsize = 22, font = :bold, rotation = π/2)

        #Collect image files
        #filter(file -> occursin(Regex("norm_$(agent_cat)\\.csv\$"),file), readdir(joinpath(pwd(),"figures","shapley",outcome)))

        images = [FileIO.load(joinpath(pwd(),"figures","shapley",outcome,img)) for img in img_files]

        fig_axes = [Axis(fig[j+1, i+1], aspect = 16/9) for i in [2,1], j in 1:3]
        #axes = [[i,j] for i in 1:2, j in 1:2]'[:]

        for (ax, img) in zip(fig_axes, images)
            image!(ax, rotr90(img))
            #axis = (aspect = DataAspect(), yreversed = true,))
            hidedecorations!(ax)
            hidespines!(ax)
        end

        # Create a legend
        colors = group_colors[1:6]
        labels = ["Flood Hazard","Flood Memory and Aversion","Location Amenities", "Housing Market", "Initial Population Distribution","Internal Stochasticity"]

        # Create legend elements for the colored items
        color_elements = [PolyElement(color = c) for c in colors]

        # Create legend element for the dashed line
        line_element = LineElement(color = :black, linestyle = :dash)

        # Combine all elements and labels
        all_elements = [color_elements..., line_element]
        all_labels = [labels..., "Flood Occurence"]

        # Add legend at the bottom (row 3, spanning both columns)
        Legend(fig[5, 2:3], all_elements, all_labels, 
            orientation = :horizontal,
            nbanks = 3,  # number of rows in the legend
            framevisible = false,
            tellwidth = false,
            tellheight = true
        )

        # Control the size of rows and columns
        rowsize!(fig.layout, 1, Auto(0.5))  # Make label row smaller
        rowsize!(fig.layout, 2, Relative(0.25))   # Make image rows larger
        rowsize!(fig.layout, 3, Relative(0.25))
        rowsize!(fig.layout, 4, Relative(0.25))
        #rowsize!(fig.layout, 4, Auto(0.5))  # Make legend row smaller

        colsize!(fig.layout, 1, Auto(0.5))  # Make label column smaller
        colsize!(fig.layout, 2, Relative(0.40))   # Make image columns larger
        colsize!(fig.layout, 3, Relative(0.40))

        #colgap!(fig.layout, 0)

        display(fig)

        CairoMakie.save(joinpath(pwd(),"figures", "shapley", "shap_indices_$(outcome)_$(scale).png"), fig)
    end
end





### Create Stacked Plot of Uncerrtainty Importance for Burdens

burd_groups = Dict(
        "Housing Market"    => ["price_inc_perc","build_inc_perc","rhea_coef","base_move"],                    
        "Location Utility"    => ["env_amen_l","env_amen_m","env_amen_h",
                                "penalty", "prop_l","prop_m","prop_h", "risk_averse", "flood_mem","flood_coefficient"
                            ],
        "Initial Population Distribution"    => ["pop_no"],
        "Internal Stochasticity"   => ["seed_order"],
        "Flood Hazard"   => ["fld_extent"], 
        "Depth-Damage" => ["DDF_order"]   
)
burd_order = ["Flood Hazard", "Housing Market", "Location Utility","Initial Population Distribution","Internal Stochasticity", "Depth-Damage"]
burd_colors = vcat(ColorSchemes.mk_8[2],ColorSchemes.mk_8[4:end])

#fld_event = "Low"
#shap_burd_ind_high = DataFrame(CSV.File(joinpath(out_dir, "post_process","shapley_indices","$(fld_event)_fld_shap_indices_flpn_burden_high.csv")))
#shap_burd_high = normalize_shap_groups(shap_burd_ind_high; normalize=true, groups = burd_groups, order = burd_order)

for outcome in ["burden", "loss"]
    #load Data and Normalize Data
    event = Int64[]
    indices = Float64[]
    stacks = Int64[]
    dodges = Int64[]

    for (fld_cat,fld_event) in enumerate(["Low", "Medium","High"])
        shap_burd_ind_high = DataFrame(CSV.File(joinpath(out_dir, "post_process","shapley_indices","$(fld_event)_fld_shap_indices_flpn_$(outcome)_high.csv")))
        shap_burd_ind_low = DataFrame(CSV.File(joinpath(out_dir, "post_process","shapley_indices","$(fld_event)_fld_shap_indices_flpn_$(outcome)_low.csv")))
        shap_burd_high = normalize_shap_groups(shap_burd_ind_high; normalize=true, groups = burd_groups, order = burd_order)
        shap_burd_low = normalize_shap_groups(shap_burd_ind_low; normalize=true, groups = burd_groups, order = burd_order)
        #Add to vectors
        append!(event,repeat([fld_cat],length(Array(shap_burd_high[1,2:end]))*2))

        append!(indices,Array(shap_burd_low[1,2:end]))
        append!(indices,Array(shap_burd_high[1,2:end]))

        append!(stacks,collect(eachindex(Array(shap_burd_low[1,2:end]))))
        append!(stacks,collect(eachindex(Array(shap_burd_high[1,2:end]))))

        append!(dodges,repeat([1],inner=length(Array(shap_burd_low[1,2:end]))),vcat(repeat([2], inner=length(Array(shap_burd_high[1,2:end])))))
    end



    fig = Figure(size = (1000, 900), fontsize = 20, pt_per_unit = 1, figure_padding = 18)
    ax = Axis(fig[1, 1], ylabel = rich("Relative Group Importance"; font = :bold), xlabel = rich("Flood Event Category"; font = :bold),
        yticks = [0.5,1.0], xgridvisible = false
    ) #
    CairoMakie.ylims!(ax, low = 0, high=1.15)
    hidespines!(ax, :t, :r)

    CairoMakie.barplot!(ax, event, indices,
            dodge = dodges,
            stack = stacks,
            color = [burd_colors[i] for i in stacks],
    )
    ax.xticks = (1:3, ["Nuisance","Moderate", "Extreme"])
    #Add labels
    CairoMakie.text!(ax, 0.65, 1.03, text=rich("Low Income", font = :italic), align = (:left, :center), rotation = (35*pi/180), fontsize = 22)
    CairoMakie.text!(ax, 1.05, 1.03, text=rich("High Income", font = :italic), align = (:left, :center), rotation = (35*pi/180), fontsize = 22)
    # Create a legend
    color_elements = [PolyElement(color = c) for c in burd_colors]

    # Add legend at the bottom (row 3, spanning both columns)
    Legend(fig[2, 1], color_elements, burd_order, 
        orientation = :horizontal,
        nbanks = 3,  # number of rows in the legend
        framevisible = false,
        tellwidth = false,
        tellheight = true
    )

    display(fig)
    CairoMakie.save(joinpath(pwd(),"figures", "shapley", "shap_indices_$(outcome).png"), fig)
end
