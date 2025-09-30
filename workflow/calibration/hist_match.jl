#activate project environment
import Pkg
Pkg.activate(".")
Pkg.instantiate()

using CSV, DataFrames
using Statistics
using ProgressMeter

include(joinpath(dirname(@__DIR__), "src", "functions.jl"))

##Set constants
threshold = 5.0
max_iterations = 10

params_df = DataFrame(CSV.File(joinpath(@__DIR__, "data/param_comb_initial_135842_ens_250.csv")))

phil_obs = DataFrame(CSV.File(joinpath(dirname(pwd()), "philadelphia-data","model_inputs", "calibration", "phil_obs_df.csv")))
obs_var = Matrix(select!(phil_obs, r"_VAR"))
###History Matching
hm_df = copy(params_df)
initial_combs = nrow(hm_df)
iteration = 1

println("Starting History Matching with $(initial_combs) parameter combinations")

println("Threshold factor (c): $threshold")
println("Max iterations: $max_iterations")
println()
mean(Matrix(select(hm_df,r"err_")), dims = 1)
mean(Matrix(select(hm_df,r"var_")), dims = 1)
#m_disc.^2 .+ ens_var.^2 .+ obs_var.^2
while iteration ≤ max_iterations
    ##Calculate average error across all combinations
    m_disc = mean(Matrix(select(hm_df,r"err_")), dims = 1)#[1,:]
    ##Calculate average variance across all combinations
    ens_var = mean(Matrix(select(hm_df,r"var_")), dims = 1)#[1,:]

    total_var = m_disc.^2 .+ ens_var.^2 .+ obs_var.^2
    ##For every combination, calculate the implausibility score for each category
    imp_cat_scores = Matrix(select(hm_df,r"err_")) ./ total_var
    imp_scores = mean(imp_cat_scores, dims=2)
    #Check which combs are below threshold
    imp_mask = imp_scores .< threshold
    #Filter param set to non-implausible combinations
    new_df = hm_df[imp_mask[:,1],:]

    #Calculate stats for this wave
    n_before = nrow(hm_df)
    n_after = nrow(new_df)
    reduction_factor = n_after / n_before
    max_implausibility = maximum(imp_scores)
    min_implausibility = minimum(imp_scores)

    #= Store wave statistics
    wave_stats = (
        wave = iteration,
        n_before = n_before,
        n_after = n_after,
        reduction_factor = reduction_factor,
        max_implausibility = max_implausibility,
        min_implausibility = min_implausibility,
        n_implausible = sum(.!non_implausible_mask)
    )
    push!(wave_statistics, wave_stats)
    =#
    println("Parameter combinations before filtering: $n_before")
    println("Parameter combinations after filtering: $n_after")
    println("Reduction factor: $(round(reduction_factor, digits=4))")
    println("Implausibility range: [$(round(min_implausibility, digits=4)), $(round(max_implausibility, digits=4))]")
    println("Implausible combinations removed: $(sum(.!imp_mask))")
    println()

    ##Check if stopping criteria is met
    if n_after == 0
        println("✓ No non-implausible parameter combinations found. Stopping.")
        break
    elseif all(imp_mask)
        println("✓ All parameter combinations are non-implausible. Stopping.")
        break
    elseif iteration > max_iterations
        println("✓ Maximum iterations ($max_iterations) reached. Stopping.")
        break
    end

    #Update for new wave/iteration
    hm_df = new_df
    println("Proceeding to next wave with $(nrow(hm_df)) combinations...")
    iteration += 1
end



# Calculate final statistics
final_stats = (
    total_waves = min(iteration, max_iterations),
    initial_combinations = initial_combs,
    final_combinations = nrow(hm_df),
    overall_reduction_factor = nrow(hm_df) / initial_combs,
    final_acceptance_rate = nrow(hm_df) / initial_combs * 100
)
    
println("=== Final Results ===")
println("Total waves completed: $(final_stats.total_waves)")
println("Initial parameter combinations: $(final_stats.initial_combinations)")
println("Final non-implausible combinations: $(final_stats.final_combinations)")
println("Overall reduction factor: $(round(final_stats.overall_reduction_factor, digits=4))")
println("Final acceptance rate: $(round(final_stats.final_acceptance_rate, digits=2))%")
println()


#Save calibrated df
CSV.write(joinpath(@__DIR__, "data/param_comb_final_135842_mean_thresh_5_ens_250.csv"), hm_df)