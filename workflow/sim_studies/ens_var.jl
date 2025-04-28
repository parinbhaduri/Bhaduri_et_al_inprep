#Determine how many ensemble members to create for each parameter combinations
#activate project environment
import Pkg
Pkg.activate(".")
Pkg.instantiate()


### PARALLEL ENSEMBLE RUN ###
include(joinpath(dirname(@__DIR__), "src", "config_parallel.jl"))
#load calibration functions
include(joinpath(dirname(@__DIR__), "src", "functions.jl"))
#Create model functions that only accept a seed value
@everywhere function phil_gen(;seed=seed)
    model = phil_model(;flood_rec = phil_flood_record, no_of_years=Int(39), start_year=Int(1981), seed=seed)      
    return model
end

gen_param = Dict(:seed => collect(range(1000,2499)))

#Run and collect data
adf,mdf = paramscan(gen_param, phil_gen; n = 39,
     parallel=true, showprogress=true, adata=simul_adata, mdata=simul_mdata
)

CSV.write(joinpath(@__DIR__,"sim_data/dataframes/adf_ens_var.csv"), adf)
CSV.write(joinpath(@__DIR__,"sim_data/dataframes/mdf_ens_var.csv"), mdf)

#Calculate model outcomes from results
using CSV, DataFrames
using DataFramesMeta
using Combinatorics

adf = DataFrame(CSV.File(joinpath(@__DIR__,"sim_data/dataframes/adf_ens_var.csv")))
mdf = DataFrame(CSV.File(joinpath(@__DIR__,"sim_data/dataframes/mdf_ens_var.csv")))

##Calculate total agent population 
adf.sum_HH = adf.sum_hh_low_HH + adf.sum_hh_med_HH + adf.sum_hh_high_HH
mdf.moved_HH = mdf.moved_low + mdf.moved_med + mdf.moved_high

##Calculate simulated outcomes of interest
#Create new Dataframe
sim_out = DataFrame(:seed => unique(mdf.seed))
#calculate total populations
pop_prop_df = @chain adf begin
    @subset(:time .== 39)
    @transform(:pop_prop_low = :sum_HH .\ :sum_hh_low_HH, :pop_prop_med = :sum_HH .\ :sum_hh_med_HH, :pop_prop_high = :sum_HH .\ :sum_hh_high_HH)
    @select(:seed, :pop_prop_low, :pop_prop_med, :pop_prop_high) 
end

sim_out = innerjoin(sim_out, pop_prop_df, on=:seed)

#Moving Proportion
move_df = @chain mdf begin
    @subset(:time .== 39)
    @transform(:moved_prop_low = :moved_HH .\ :moved_low, :moved_prop_med = :moved_HH .\ :moved_med, :moved_prop_high = :moved_HH .\ :moved_high)
    @select(:seed, :moved_prop_low, :moved_prop_med, :moved_prop_high) 
end

sim_out = innerjoin(sim_out, move_df, on=:seed)

#Pop Growth and Sales Price Growth
pop_price_change_df = @chain adf begin
    @subset((:time .== 10) .| (:time .== 39))
    @groupby(:seed)
    @combine(:pop_growth_low = (:sum_hh_low_HH[2] - :sum_hh_low_HH[1])/:sum_hh_low_HH[1], :pop_growth_med = (:sum_hh_med_HH[2] - :sum_hh_med_HH[1])/:sum_hh_med_HH[1],
     :pop_growth_high = (:sum_hh_med_HH[2] - :sum_hh_high_HH[1])/:sum_hh_high_HH[1], :price_growth_low = (:mean_price_low_BG[2] - :mean_price_low_BG[1])/:mean_price_low_BG[1],
     :price_growth_med = (:mean_price_med_BG[2] - :mean_price_med_BG[1])/:mean_price_med_BG[1],
     :price_growth_high = (:mean_price_high_BG[2] - :mean_price_high_BG[1])/:mean_price_high_BG[1])
    #@transform(:moved_prop_low = :moved_HH .\ :moved_low, :moved_prop_med = :moved_HH .\ :moved_med, :moved_prop_high = :moved_HH .\ :moved_high)
    #@select(:seed, :moved_prop_low, :moved_prop_med, :moved_prop_high) 
end


sim_out = innerjoin(sim_out, pop_price_change_df, on=:seed)
sim_out_low = select(sim_out, r"_low")
sim_out_med = select(sim_out, r"_med")
sim_out_high = select(sim_out, r"_high")

## Calculate Ensemble Variance of model errors for different ensemble sizes
ens_size = collect(range(10, 1500, step=10))
ens_var_low = []
ens_var_med = []
ens_var_high = []
for e_s in ens_size
    push!(ens_var_low, var([model_error(Vector(sim_out_low[iter[1],:]),Vector(sim_out_low[iter[2],:])) for iter in Combinatorics.combinations(round.(Int, range(1,e_s, step = 1)), 2)])) #round.(Int, range(1,1500, length = e_s))
    push!(ens_var_med, var([model_error(Vector(sim_out_med[iter[1],:]),Vector(sim_out_med[iter[2],:])) for iter in Combinatorics.combinations(round.(Int, range(1,e_s, step = 1)), 2)]))
    push!(ens_var_high, var([model_error(Vector(sim_out_high[iter[1],:]),Vector(sim_out_high[iter[2],:])) for iter in Combinatorics.combinations(round.(Int, range(1,e_s, step = 1)), 2)]))
end

##Plot results
using Plots
using ColorSchemes
include("sim_functions.jl")
### Plot Total ensemble
ens_plots = simul_plot(adf, :seed; leg = false, color = cgrad(:grays))
ens_mark_plots = simul_market(adf,mdf, :seed; leg = false, color = cgrad(:grays))
#Add averages
gadf = groupby(adf, :time)
avg_adf = combine(gadf, Not(:time) .=> mean, renamecols=false)
gmdf = groupby(mdf, :time)
avg_mdf = combine(gmdf, Not(:time) .=> mean, renamecols=false)

for (i,col) in enumerate([:sum_hh_low_HH,:sum_occ_low_BG,:sum_hh_med_HH,:sum_occ_med_BG,:sum_hh_high_HH,:sum_occ_high_BG])
    plot!(ens_plots[i], avg_adf.time, avg_adf[!,col], lw = 3, lc = :red)
end

for (i,col) in enumerate([:moved_low,:mean_price_low_BG,:moved_med,:mean_price_med_BG,:moved_high,:mean_price_high_BG])
    if isodd(i)
        plot!(ens_mark_plots[i], avg_mdf.time[2:end], avg_mdf[!,col][2:end], lw = 3, lc = :red)
    else
        plot!(ens_mark_plots[i], avg_adf.time, avg_adf[!,col], lw = 3, lc = :red)
    end
end


plot(ens_plots..., layout=(3, 2), size = (1100, 1000), plot_title = "Pop. Dynamics (Ensemble)")
plot(ens_mark_plots..., layout=(3, 2), size = (1100, 1000), plot_title = "Market Dynamics (Ensemble)")

### Plot changes in ensemble variance
plot(ens_size, [ens_var_low ens_var_med ens_var_high], label = ["Low Income" "Middle Income" "High Income"], lw=3)
xlabel!("Ensemble Size")
ylabel!("Variance in Error")