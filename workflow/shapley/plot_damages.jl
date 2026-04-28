# load environment and packages
import Pkg
Pkg.activate(".")
Pkg.instantiate()

using CSV, DataFrames # read CSV of Shapley indices
using CairoMakie # plotting library
using ColorSchemes
using Measures # adjust margins with explicit measures
using Statistics # get mean function
using StatsPlots
using HDF5
using Parquet2
#using CategoricalArrays


### Load Data ###
out_dir = joinpath(@__DIR__,"data","shap_runs")

Palette = ColorSchemes.Hokusai3
hc = Palette[2] #High Income Color
mc = Palette[1] #Middle Income Color
lc = Palette[4] #low Income Color

#Initialize 2 plots: One for total estimates and One for total discrepancies
fig1 = Figure(size = (1000, 1000), fontsize = 24, pt_per_unit = 1, figure_padding = 18)
fig2 = Figure(size = (1000, 1000), fontsize = 24, pt_per_unit = 1, figure_padding = 18)
ga = fig1[1, 1:2] = GridLayout()
gb = fig1[2, 1:2] = GridLayout()

gc = fig2[1, 1:2] = GridLayout()
gd = fig2[2, 1:2] = GridLayout()

ax1_events = Axis(ga[1,1],
    xaxisposition = :top,
    #xticklabelpad = 50.0,
    xticks = ([2,5,8], string.(["Nuisance","Moderate", "Extreme"])), xticksvisible = false,
    xlabel = rich("Flood Event Category"; font = :bold), xgridvisible = false
)

hidespines!(ax1_events)
hideydecorations!(ax1_events)

ax1 = Axis(ga[1, 1], ylabel = rich("Flood Loss (Million USD)"; font = :bold), xlabel = "",#rich("Year of Event"; font = :bold),
    yticks = ([0.0,0.5e8,1.0e8,1.5e8], string.([0,50,100,150])),  
    xticks = (collect(1:9), string.(collect(1:9))),#["2010","2013","1988","1981","2018","1991","1996","1989","2011"]),
    xticklabelsize = 20, xlabelsize = 22, xgridvisible = false, xticksvisible = true#, xticklabelsvisible = false,
) #
hidespines!(ax1, :t, :r)

ax2 = Axis(gb[1, 1], ylabel = rich("Flood Loss Burden (%)"; font = :bold), xlabel = rich("Flood Event"; font = :bold),
    #yticks = ([1.0e8,2.0e8,3.0e8], string.([100,200,300])),
    xticks = (collect(1:9), string.(collect(1:9))),#["2010","2013","1988","1981","2018","1991","1996","1989","2011"]),
    xticklabelsize = 20, xlabelsize = 22, xgridvisible = false, #xticklabelsvisible = false, xticksvisible = false
) #
#ylims!(ax2, low = 0, high=20)
hidespines!(ax2, :t, :r)


ax3_events = Axis(gc[1,1],
    xaxisposition = :top,
    #xticklabelpad = 50.0,
    xticks = ([2,5,8], string.(["Nuisance","Moderate", "Extreme"])), xticksvisible = false,
    xlabel = rich("Flood Event Category"; font = :bold), xgridvisible = false
)

hidespines!(ax3_events)
hideydecorations!(ax3_events)

ax3 = Axis(gc[1, 1], ylabel = rich("Total Flood Loss\nDiscrepancy (Million USD)"; font = :bold), xlabel = "",#rich("Year of Event"; font = :bold),
    yticks = ([0.0,1.0e7,2.0e7,3.0e7], string.([0,10,20,30])), 
    xticks = (collect(1:9), string.(collect(1:9))),#["2010","2013","1988","1981","2018","1991","1996","1989","2011"]),
    xticklabelsize = 20, xlabelsize = 22, xgridvisible = false, xticksvisible = true#, xticklabelsvisible = false,
) #
hidespines!(ax3, :t, :r)

ax4 = Axis(gd[1, 1], ylabel = rich("Total Flood Loss Burden\nDiscrepancy (%)"; font = :bold), xlabel = rich("Flood Event"; font = :bold),
    yticks = ([-0.20,-0.10,0,0.10,0.20], string.([-20,-10,0,10,20])), 
    xticks = (collect(1:9), string.(collect(1:9))),#["2010","2013","1988","1981","2018","1991","1996","1989","2011"]),
    xticklabelsize = 20, xlabelsize = 22, xgridvisible = false, xticksvisible = true#, xticklabelsvisible = false,
) #
hidespines!(ax4, :t, :r)




data_len = 159500
groups = repeat([1,2,3],inner=data_len*3)
dodges = repeat(vcat(repeat([1], inner=data_len),repeat([2],inner=data_len),repeat([3],inner=data_len)),outer=3)
colors = vcat(
    fill(lc, data_len),
    fill(mc, data_len),
    fill(hc, data_len),
    fill(lc, data_len),
    fill(mc, data_len),
    fill(hc, data_len),
    fill(lc, data_len),
    fill(mc, data_len),
    fill(hc, data_len)   
);
# Create dict of Haz categories and flood years
Haz_Dict = Dict(
    "High" => [1996,1989,2011],
    "Medium" => [1981,2018,1991],
    "Low" => [2010,2013,1988],
)

for (ind,haz_size) in enumerate(["Low", "Medium", "High"])
    loss_data = []
    loss_diff = []
    burden_data = []
    burd_diff = []
    sampled_indices_l = range(1, 79909500, length=159500) |> x -> round.(Int, x) #includes default assumptions (no uncertainty) distribution
    sampled_indices_s = range(1, 79750000, length=159500) |> x -> round.(Int, x)
    for year in Haz_Dict[haz_size]
        println("Starting year $(year)...")
        #First, load damages with no uncertainty
        dam_file = h5open(joinpath(out_dir, "post_process","flood_loss",haz_size,"$(year)_flood_loss.h5"), "r")
        dam_no_unc = dam_file["no_unc estimates"]
        flpn_loss_no_unc = zeros(size(dam_no_unc,1),3,2)
        for inc in [4,5,6] 
            for i in 1:1000:size(dam_no_unc,1)
                end_idx = min(i + 1000 - 1, size(dam_no_unc,1))
                chk = dam_no_unc[i:end_idx,:,inc]
                price_chk = dam_no_unc[i:end_idx,:,inc-3]
                loss_sum = dropdims(sum(x -> isnan(x) ? 0 : x, chk, dims=2), dims=2)
                #Calculate burden
                tot_price = dropdims(sum(x -> isnan(x) ? 0 : x, price_chk, dims=2), dims=2)
                outcome = loss_sum ./ tot_price
                replace!(outcome, NaN => 0) #NaN values occur when no agents exist in area. Assume burden is 0
                flpn_loss_no_unc[i:end_idx,inc-3,1] = loss_sum
                flpn_loss_no_unc[i:end_idx,inc-3,2] = outcome
            end
        end
        #Convert no_unc loss to df
        no_unc_df = DataFrame(:no_unc_loss_low => flpn_loss_no_unc[:,1,1],:no_unc_loss_med => flpn_loss_no_unc[:,2,1], :no_unc_loss_high => flpn_loss_no_unc[:,3,1],
                            :no_unc_burd_low => flpn_loss_no_unc[:,1,2],:no_unc_burd_med => flpn_loss_no_unc[:,2,2], :no_unc_burd_high => flpn_loss_no_unc[:,3,2],
                            :model => 1:size(flpn_loss_no_unc,1)
        )

        #Load loss values with uncertainty
        ds = Parquet2.Dataset(joinpath(out_dir, "post_process","flood_loss",haz_size))
        filtered_file_high = findfirst(file -> occursin(Regex("$(year)_model_outcome_flpn_loss_high"),file), string.(Parquet2.filelist(ds)))
        append!(ds,filtered_file_high)
        filtered_file_med = findfirst(file -> occursin(Regex("$(year)_model_outcome_flpn_loss_med"),file), string.(Parquet2.filelist(ds)))
        append!(ds,filtered_file_med)
        filtered_file_low = findfirst(file -> occursin(Regex("$(year)_model_outcome_flpn_loss_low"),file), string.(Parquet2.filelist(ds)))
        append!(ds,filtered_file_low)
        filtered_file_total = findfirst(file -> occursin(Regex("$(year)_model_outcome_flpn_loss.parquet"),file), string.(Parquet2.filelist(ds)))
        append!(ds,filtered_file_total)

        event_loss_inc_high = ds[1] |> Parquet2.select(:model, :loss_value, :burden_value) |> DataFrame
        event_loss_inc_high = innerjoin(event_loss_inc_high, no_unc_df[:,[:no_unc_loss_high, :no_unc_burd_high, :model]], on = :model)
        event_loss_high = sort(vcat(no_unc_df.no_unc_loss_high, event_loss_inc_high.loss_value))
        event_burden_high = sort(vcat(no_unc_df.no_unc_burd_high, event_loss_inc_high.burden_value))
        event_diff_high = sort(event_loss_inc_high.no_unc_loss_high .- event_loss_inc_high.loss_value)
        event_diff_burd_high = sort(event_loss_inc_high.no_unc_burd_high .- event_loss_inc_high.burden_value)
        
        event_loss_inc_med = ds[2] |> Parquet2.select(:model, :loss_value, :burden_value) |> DataFrame
        event_loss_inc_med = innerjoin(event_loss_inc_med, no_unc_df[:,[:no_unc_loss_med, :no_unc_burd_med, :model]], on = :model)
        event_loss_med = sort(vcat(no_unc_df.no_unc_loss_med, event_loss_inc_med.loss_value))
        event_burden_med = sort(vcat(no_unc_df.no_unc_burd_med, event_loss_inc_med.burden_value))
        event_diff_med = sort(event_loss_inc_med.no_unc_loss_med .- event_loss_inc_med.loss_value)
        event_diff_burd_med = sort(event_loss_inc_med.no_unc_burd_med .- event_loss_inc_med.burden_value)
        
        event_loss_low = Parquet2.load(ds[3],"loss_value")
        event_loss_inc_low = ds[3] |> Parquet2.select(:model, :loss_value, :burden_value) |> DataFrame
        event_loss_inc_low = innerjoin(event_loss_inc_low, no_unc_df[:,[:no_unc_loss_low, :no_unc_burd_low, :model]], on = :model)
        event_loss_low = sort(vcat(no_unc_df.no_unc_loss_low, event_loss_inc_low.loss_value))
        event_burden_low = sort(vcat(no_unc_df.no_unc_burd_low, event_loss_inc_low.burden_value))
        event_diff_low = sort(event_loss_inc_low.no_unc_loss_low .- event_loss_inc_low.loss_value)
        event_diff_burd_low = sort(event_loss_inc_low.no_unc_burd_low .- event_loss_inc_low.burden_value)
        
        
        #Add sets to data vectors
        append!(loss_data,vcat(event_loss_low[sampled_indices_l], event_loss_med[sampled_indices_l], event_loss_high[sampled_indices_l]))
        append!(burden_data,vcat(event_burden_low[sampled_indices_l].*100, event_burden_med[sampled_indices_l].*100, event_burden_high[sampled_indices_l].*100))
        append!(loss_diff,vcat(event_diff_low[sampled_indices_s], event_diff_med[sampled_indices_s], event_diff_high[sampled_indices_s]))
        append!(burd_diff,vcat(event_diff_burd_low[sampled_indices_s], event_diff_burd_med[sampled_indices_s], event_diff_burd_high[sampled_indices_s]))
        
    end
    println("Plotting Boxplots for Hazard Category $(haz_size)")
    CairoMakie.boxplot!(ax1, groups.+((ind-1)*3), loss_data, 
        dodge = dodges,
        color = colors,
        strokewidth = 1,
        #width = 0.5,  # Make boxes narrower
    )

    CairoMakie.boxplot!(ax2, groups.+((ind-1)*3), burden_data, 
        dodge = dodges,
        color = colors,
        strokewidth = 1,
        #width = 0.5,  # Make boxes narrower
    )

    CairoMakie.boxplot!(ax3, groups.+((ind-1)*3), loss_diff, 
        dodge = dodges,
        color = colors,
        strokewidth = 1,
        #width = 0.5,  # Make boxes narrower
    )

    CairoMakie.boxplot!(ax4, groups.+((ind-1)*3), burd_diff, 
        dodge = dodges,
        color = colors,
        strokewidth = 1,
        #width = 0.5,  # Make boxes narrower
    )
    
end
CairoMakie.vlines!(ax1, [3.5,6.5], color = :grey, linecap = :round, linewidth = 2.5)  #linestyle = :dash, ymax = 0.8,
CairoMakie.vlines!(ax2, [3.5,6.5], color = :grey, linecap = :round, linewidth = 2.5) 
CairoMakie.vlines!(ax3, [3.5,6.5], color = :grey, linecap = :round, linewidth = 2.5)
CairoMakie.vlines!(ax4, [3.5,6.5], color = :grey, linecap = :round, linewidth = 2.5)
#Add arrow for extent direction
text!(ax1, 3.65, 1.6e8, text=rich("smallest to largest flood extent ", font = :italic), align = (:left, :center), fontsize = 18)
arrows!(ax1, [4.0], [1.5e8], [2.0], [0.0], linewidth=2.0)
text!(ax3, 3.65, 3.2e7, text=rich("smallest to largest flood extent ", font = :italic), align = (:left, :center), fontsize = 18)
arrows!(ax3, [4.0], [3.0e7], [2.0], [0.0], linewidth=2.0)
#Create Legend
elem_1 = [PolyElement(color = lc)]
elem_2 = [PolyElement(color = mc)]
elem_3 = [PolyElement(color = hc)]

axislegend(ax1, [elem_1, elem_2, elem_3] , ["Low Income", "Middle Income", "High Income"],
    orientation = :vertical, framevisible = false, labelsize = 22, position = :lt
)

axislegend(ax3, [elem_1, elem_2, elem_3] , ["Low Income", "Middle Income", "High Income"],
    orientation = :vertical, framevisible = false, labelsize = 22, position = :lt
)

rowgap!(fig1.layout, 50)
rowgap!(fig2.layout, 50)
display(fig1)
CairoMakie.save(joinpath(pwd(),"figures", "loss", "flpn_flood_loss_abs.png"), fig1)
display(fig2)
CairoMakie.save(joinpath(pwd(),"figures", "loss", "flpn_flood_loss_discrep.png"), fig2)




## Plot for Damage Uncertainty
tc = ColorSchemes.Archambault[2] #Total Damages color
tc_no_unc = ColorSchemes.Archambault[3]

fig = Figure(size = (1000, 800), fontsize = 24, pt_per_unit = 1, figure_padding = 18)
ga = fig[1, 1:2] = GridLayout()


ax1_events = Axis(ga[1,1],
    xaxisposition = :top,
    #xticklabelpad = 50.0,
    xticks = ([2,5,8], string.(["0-50%","51-90%", "91-100%"])), xticksvisible = false,
    xlabel = rich("Flood Hazard Extent"; font = :bold), xgridvisible = false
)
hidespines!(ax1_events)
hideydecorations!(ax1_events)

ax1 = Axis(ga[1, 1], ylabel = rich("Flood Loss\n(Million USD)"; font = :bold), xlabel = rich("Year of Event"; font = :bold),
    yticks = ([0.00,1.0e8,2.0e8,3.0e8], string.([0,100,200,300])), #xticks = (collect(1:9), repeat(["High", "Middle","Low"],3)),
    xticks = (collect(1:9), ["2010","2013","1988","1981","2018","1991","1996","1989","2011"]),
    xticklabelsize = 20, xlabelsize = 22, xgridvisible = false, #xticklabelsvisible = false, xticksvisible = false
) #
hidespines!(ax1, :t, :r)



#=
#Check raw damages
phil_damages = DataFrame(CSV.File(joinpath(dirname(pwd()), "philadelphia-data","flood_hazard", "data","phil_flood_dmg_ens.csv")))
phil_damages[phil_damages.bg_id .== 421010087013, :naccs_loss_2011]
combine(groupby(phil_damages,:bg_id), All() .=> sum)
phil_damages_no_unc = DataFrame(CSV.File(joinpath(dirname(pwd()), "philadelphia-data","flood_hazard", "data","phil_flood_dmg_no_unc.csv")))

no_unc_tot = combine(groupby(phil_damages_no_unc, :bg_id), All() .=> sum)
no_unc_tot[no_unc_tot.naccs_loss_2011_sum .== maximum(no_unc_tot.naccs_loss_2011_sum),:]
=#

unc_groups = repeat([1,2,3],inner=data_len*2)
unc_dodges = repeat(vcat(repeat([1], inner=data_len),repeat([2],inner=data_len)),outer=3)
unc_colors = vcat(
    fill(tc, data_len),
    fill(tc_no_unc, data_len),
    fill(tc, data_len),
    fill(tc_no_unc, data_len),
    fill(tc, data_len),
    fill(tc_no_unc, data_len)   
);


for (ind,haz_size) in enumerate(["Low", "Medium", "High"])
    total_loss = []
    sampled_indices = range(1, 79750000, length=159500) |> x -> round.(Int, x)
    for year in Haz_Dict[haz_size]
        println("Starting year $(year)...")
        #First, load damages with no uncertainty
        dam_file = h5open(joinpath(out_dir, "post_process","flood_loss",haz_size,"$(year)_flood_loss.h5"), "r")
        dam_no_unc = dam_file["no_unc estimates"]
        flpn_loss_no_unc = zeros(size(dam_no_unc,1),3)
        for inc in [4,5,6] 
            for i in 1:1000:size(dam_no_unc,1)
                end_idx = min(i + 1000 - 1, size(dam_no_unc,1))
                chk = dam_no_unc[i:end_idx,:,inc]
                loss_sum = dropdims(sum(x -> isnan(x) ? 0 : x, chk, dims=2), dims=2)
                flpn_loss_no_unc[i:end_idx,inc-3] .= loss_sum
            end
        end
        #Convert no_unc loss to df
        no_unc_df = DataFrame(:no_unc_value_low => flpn_loss_no_unc[:,1],:no_unc_value_med => flpn_loss_no_unc[:,2], :no_unc_value_high => flpn_loss_no_unc[:,3],:model => 1:length(flpn_loss_no_unc[:,1]))

        #Load loss values with uncertainty
        ds = Parquet2.Dataset(joinpath(out_dir, "post_process","flood_loss",haz_size))
        filtered_file_total = findfirst(file -> occursin(Regex("$(year)_model_outcome_flpn_loss.parquet"),file), string.(Parquet2.filelist(ds)))
        append!(ds,filtered_file_total)
        event_loss_total = Parquet2.load(ds[1],"loss_value")
        sort!(event_loss_total)
        
        #Add both sets to data vectors
        append!(total_loss,vcat(event_loss_total[sampled_indices], dropdims(sum(x -> isnan(x) ? 0 : x, flpn_loss_no_unc, dims=2), dims=2)))
    end
    println("Plotting Boxplots for Hazard Category $(haz_size)")    
    CairoMakie.boxplot!(ax1, unc_groups.+((ind-1)*3), total_loss, 
        dodge = unc_dodges,
        color = unc_colors,
        strokewidth = 1,
        #width = 0.5,  # Make boxes narrower
    )
    
end

CairoMakie.vlines!(ax1, [3.5,6.5], color = :grey, linecap = :round, linewidth = 2.5)  #linestyle = :dash, ymax = 0.8,
 
#Add arrow for extent direction
text!(ax1, 3.65, 2.4e8, text=rich("largest to smallest flood extent ", font = :italic), align = (:left, :center), fontsize = 18)
arrows!(ax1, [4.0], [2.2e8], [2.0], [0.0], linewidth=2.0)

#Create Legend for total losses
elem_1 = [PolyElement(color = tc)]
elem_2 = [PolyElement(color = tc_no_unc)]

axislegend(ax1, [elem_1, elem_2] , ["Damage Uncertainty", "Default Assumptions"],
    orientation = :vertical, framevisible = false, labelsize = 22, position = :lt
)

display(fig)
CairoMakie.save(joinpath(pwd(),"figures", "loss", "flpn_flood_loss_unc.png"), fig)




