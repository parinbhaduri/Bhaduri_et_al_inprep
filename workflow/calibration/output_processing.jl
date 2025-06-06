#activate project environment
import Pkg
Pkg.activate(".")
Pkg.instantiate()

using CSV, DataFrames
using DataFramesMeta

#Load data
#t_d = DataFrame(CSV.File(joinpath(@__DIR__,"data/results_118418/agents_chunk_2.csv")))
dat_dir =joinpath(@__DIR__,"data/results_118418/")

param_cols = ["env_amen_l","price_inc_perc","rhea_coef","base_move",
"prop_l","prop_m","env_amen_m","risk_averse","build_inc_perc",
"env_amen_h","flood_coefficient","penalty","prop_h", "seed"]


agent_files = filter(file -> occursin(r"^agents.*\.csv$",file), readdir(dat_dir))
model_files = filter(file -> occursin(r"^model.*\.csv$",file), readdir(dat_dir))


model_params_df = DataFrame([name => Float64[] for name in param_cols])
agent_params_df = DataFrame([name => Float64[] for name in param_cols])

for file in model_files
    df = DataFrame(CSV.File(joinpath(dat_dir,file)))
    df.moved_HH = df.moved_low + df.moved_med + df.moved_high
    #Moving Proportion
    move_df = @chain df begin
        @subset(:time .== 39)
        @transform(:moved_prop_low = :moved_HH .\ :moved_low, :moved_prop_med = :moved_HH .\ :moved_med, :moved_prop_high = :moved_HH .\ :moved_high)
        @select(Between(:env_amen_l, :moved_prop_high)) 
    end
    append!(model_params_df, move_df, cols = :union)
    move_df = nothing
    df = nothing
end

params_seeds_df

for file in agent_files
    df = DataFrame(CSV.File(joinpath(dat_dir,file)))
    df.sum_HH = df.sum_hh_low_HH + df.sum_hh_med_HH + df.sum_hh_high_HH
    #calculate total populations
    pop_prop_df = @chain df begin
        @subset(:time .== 39)
        @transform(:pop_prop_low = :sum_HH .\ :sum_hh_low_HH, :pop_prop_med = :sum_HH .\ :sum_hh_med_HH, :pop_prop_high = :sum_HH .\ :sum_hh_high_HH)
        @select(Between(:env_amen_l, :pop_prop_high)) 
    end
    #Pop Growth and Sales Price Growth
    pop_price_change_df = @chain df begin
        @subset((:time .== 20) .| (:time .== 39))
        @groupby(param_cols)
        @combine(:pop_growth_low = (:sum_hh_low_HH[2] - :sum_hh_low_HH[1])/:sum_hh_low_HH[1], :pop_growth_med = (:sum_hh_med_HH[2] - :sum_hh_med_HH[1])/:sum_hh_med_HH[1],
        :pop_growth_high = (:sum_hh_med_HH[2] - :sum_hh_high_HH[1])/:sum_hh_high_HH[1], :price_growth_low = (:mean_price_low_BG[2] - :mean_price_low_BG[1])/:mean_price_low_BG[1],
        :price_growth_med = (:mean_price_med_BG[2] - :mean_price_med_BG[1])/:mean_price_med_BG[1],
        :price_growth_high = (:mean_price_high_BG[2] - :mean_price_high_BG[1])/:mean_price_high_BG[1])
        #@select(:seed, :moved_prop_low, :moved_prop_med, :moved_prop_high) 
    end
    #Join both agent attribute dfs together
    adf = innerjoin(pop_prop_df, pop_price_change_df, on=param_cols)
    #Join to main df
    append!(agent_params_df, adf, cols = :union)

    pop_prop_df = nothing
    pop_price_change_df = nothing
    adf = nothing
    df = nothing
end

params_df = innerjoin(model_params_df, agent_params_df, on=param_cols)

CSV.write(joinpath(@__DIR__, "data/calib_sim_output.csv"), params_df)









###test###
t_df = DataFrame(CSV.File(joinpath(dat_dir,model_files[1])))
#Calculate total moving population for all time steps
t_df.moved_HH = t_df.moved_low + t_df.moved_med + t_df.moved_high
t_move_df = @chain t_df begin
    @subset(:time .== 39)
    @transform(:moved_prop_low = :moved_HH .\ :moved_low, :moved_prop_med = :moved_HH .\ :moved_med, :moved_prop_high = :moved_HH .\ :moved_high)
    @select(Between(:env_amen_l, :moved_prop_high))
    #@select(vcat(Symbol.(param_cols), :moved_prop_low, :moved_prop_med, :moved_prop_high))
end

#calculate total populations
t_a_df = DataFrame(CSV.File(joinpath(dat_dir,agent_files[1])))
t_a_df.sum_HH = t_a_df.sum_hh_low_HH + t_a_df.sum_hh_med_HH + t_a_df.sum_hh_high_HH
pop_prop_df = @chain t_a_df begin
    @subset(:time .== 39)
    @transform(:pop_prop_low = :sum_HH .\ :sum_hh_low_HH, :pop_prop_med = :sum_HH .\ :sum_hh_med_HH, :pop_prop_high = :sum_HH .\ :sum_hh_high_HH)
    @select(Between(:env_amen_l, :pop_prop_high)) 
end

#Pop Growth and Sales Price Growth
pop_price_change_df = @chain t_a_df begin
    @subset((:time .== 10) .| (:time .== 39))
    @groupby(param_cols)
    @combine(:pop_growth_low = (:sum_hh_low_HH[2] - :sum_hh_low_HH[1])/:sum_hh_low_HH[1], :pop_growth_med = (:sum_hh_med_HH[2] - :sum_hh_med_HH[1])/:sum_hh_med_HH[1],
     :pop_growth_high = (:sum_hh_med_HH[2] - :sum_hh_high_HH[1])/:sum_hh_high_HH[1], :price_growth_low = (:mean_price_low_BG[2] - :mean_price_low_BG[1])/:mean_price_low_BG[1],
     :price_growth_med = (:mean_price_med_BG[2] - :mean_price_med_BG[1])/:mean_price_med_BG[1],
     :price_growth_high = (:mean_price_high_BG[2] - :mean_price_high_BG[1])/:mean_price_high_BG[1])
    #@transform(:moved_prop_low = :moved_HH .\ :moved_low, :moved_prop_med = :moved_HH .\ :moved_med, :moved_prop_high = :moved_HH .\ :moved_high)
    #@select(:seed, :moved_prop_low, :moved_prop_med, :moved_prop_high) 
end

innerjoin(pop_prop_df, pop_price_change_df, on=param_cols)