### Calculate output metrics from model run data ###

# load environment and packages
import Pkg
Pkg.activate(".")
Pkg.instantiate()


using CSV, DataFrames
using Statistics 
using HDF5
using ProgressMeter

###Read in data###
out_dir = joinpath(@__DIR__,"data","shap_runs")
filename = "2011_abm_data_142772.h5"

h5file = h5open(joinpath(out_dir,filename), "r")
pop_dat = h5file["pop_data"]
price_dat = h5file["price_data"]
println("Pop Dataset size: ", size(pop_dat))
println(h5file["pop_vars"][:])
println("Price Dataset size: ", size(price_dat))
println(h5file["price_vars"][:])
#println("Data type: ", eltype(dataset))

#Subset Area to Floodplain (Select Block Groups)
read(h5file["historical flood year"])
exp_bgs = phil_flood_record[phil_flood_record[!,string(read(h5file["historical flood year"]))] .> 0,"GEOID"]
flpn_index = filter(x -> x !== nothing, indexin(exp_bgs, h5file["GEOID"][:]))


chunk_size = 15950  # Adjust based on your RAM
total_rows = Int(size(pop_dat, 1))
        
println("Processing $total_rows rows in chunks of $chunk_size")
function process_output(dataset;chunk_size=6050,total_rows=159500, var_col = 1, subset=false, index=flpn_index)
    output = zeros(total_rows,size(dataset, 3))
    for i in 1:chunk_size:total_rows
        end_idx = min(i + chunk_size - 1, total_rows)
        chunk = dataset[i:end_idx, :, :, var_col]
        if subset
            outcome = zeros(size(chunk, 1),size(chunk, 3))
            for ind in index
                outcome += chunk[:,ind,:]
            end
        else
            outcome = dropdims(sum(chunk,dims=2), dims=2)
        end
        
        # Add chunk data to output
        output[i:end_idx, :] = outcome
        # Clear memory
        flpn_low = nothing
        chunk = nothing
        GC.gc()
    end 
    return output
end



flpn_low = process_output(pop_dat;var_col=4, subset=true)
flpn_med = process_output(pop_dat;var_col=5, subset=true)
flpn_high = process_output(pop_dat;var_col=6, subset=true)

close(h5file)


init_pop_lo = flpn_low[:, 1]
flpn_norm_low = (flpn_low .- init_pop_lo) ./ init_pop_lo .* 100

init_pop_hi = flpn_high[:, 1]
flpn_norm_high = (flpn_high .- init_pop_hi) ./ init_pop_hi .* 100

#Repeat for entire floodplain pop 
pop_flpn = flpn_low .+ flpn_med .+ flpn_high
init_pop = pop_flpn[:, 1]
flpn_norm = (pop_flpn .- init_pop) ./ init_pop .* 100

#pop_flpn = flpn_low .+ flpn_med .+ flpn_high
#low_prop = flpn_low ./ pop_flpn

out_df = DataFrame(flpn_norm,Symbol.(1979 .+ collect(1:size(flpn_norm,2))))
CSV.write(joinpath(out_dir, "post_process","population","model_outcome_flpn_pop_norm.csv"), out_df)

#Calculate median and 95% intervals
time = collect(range(0.0,39.0))
low_med = mapslices(x -> median(skipmissing(x)), flpn_norm_low, dims=1)
low_quantiles = mapslices(x -> quantile(skipmissing(x), [0.025, 0.975]), flpn_norm_low, dims=1)

hi_med = mapslices(x -> median(skipmissing(x)), flpn_norm_high, dims=1)
hi_quantiles = mapslices(x -> quantile(skipmissing(x), [0.025, 0.975]), flpn_norm_high, dims=1)

#Plot results
using CairoMakie
using ColorSchemes

Palette = ColorSchemes.Hokusai3

fig = Figure(size = (1000,800), fontsize = 24, pt_per_unit = 1, figure_padding = 18)

ax1 = Axis(fig[1, 1], ylabel = rich("Change in Population (%)";font=:bold), xlabel = rich("Time (years)";font=:bold), #xscale = log10,
titlesize = 28,  xgridvisible = false, titlealign = :center, title = "Change in Population in Exposed Areas") #xticks = ([10,100,1000], string.([10,100,1000])), limits = ((10,1000), nothing),

CairoMakie.lines!(ax1, time, vec(low_med), color = Palette[4], linewidth = 3.5)
#, label = false)
CairoMakie.band!(ax1, time, low_quantiles[1,:], low_quantiles[2,:], color = (Palette[4], 0.35))

#display(fig)


#fig = Figure(size = (1000,800), fontsize = 16, pt_per_unit = 1, figure_padding = 18)

#ax1 = Axis(fig[1, 1], ylabel = rich("Change in Population (%)";font=:bold), xlabel = rich("Time (years)";font=:bold), #xscale = log10,
#titlesize = 18,  xgridvisible = false, titlealign = :left, title = "High Income Population in Floodplain") #xticks = ([10,100,1000], string.([10,100,1000])), limits = ((10,1000), nothing),

CairoMakie.lines!(ax1, time, vec(hi_med), color = Palette[2], linewidth = 3.5)
#, label = false)
CairoMakie.band!(ax1, time, hi_quantiles[1,:], hi_quantiles[2,:], color = (Palette[2], 0.35))


#Create Legend
elem_1 = [LineElement(color = Palette[1], linestyle = :solid, linewidth = 5), PolyElement(color = (Palette[4], 0.35))]

elem_2 = [LineElement(color = Palette[2], linestyle = :solid, linewidth = 5), PolyElement(color = (Palette[2], 0.35))]

axislegend(ax1, [elem_1, elem_2] , ["Low Income Population", "High Income Population"],
orientation = :vertical, framevisible = false, labelsize = 26, position = :lt) #

display(fig)

CairoMakie.save(joinpath(pwd(),"figures", "shapley", "2011_pop_recovery.png"), fig)