#activate project environment
import Pkg
Pkg.activate(".")
Pkg.instantiate()


using Distributed, SlurmClusterManager


addprocs(SlurmManager())

#using Distributed
#addprocs(12, exeflags="--project=$(Base.active_project())")

# instantiate and precompile environment
@everywhere begin
  using Pkg;Pkg.activate(".");
  Pkg.instantiate(); Pkg.precompile()
end


### PARALLEL ENSEMBLE RUN ###
@everywhere begin
    using Dates
    using ProgressMeter
    using CSV, DataFrames
    using Statistics
    using DataStructures
    using Agents
    using CHANCE_C
    using LinearAlgebra
    using HDF5
end


@everywhere include(joinpath(dirname(@__DIR__),"src/data_collect.jl"))
@everywhere include("shap_functions.jl")




#Load calibrated parameter combinations
param_path = joinpath(dirname(@__DIR__),"calibration","data/param_comb_final_mean_thresh_6_ens_250.csv")
calib_combs = DataFrame(CSV.File(param_path))[:,1:14]

#Load flood hazard categories
haz_cat = DataFrame(CSV.File(joinpath(dirname(pwd()), "philadelphia-data","model_inputs", "phil_flood_hist_categories.csv")))

events = combine(groupby(haz_cat, "category")) do group
    # Find indices for min and max flood extents
    min_idx = argmin(group.total_extents)
    max_idx = argmax(group.total_extents)
    
    # For median, sort and find middle index
    sorted_indices = sortperm(group.total_extents)
    median_idx = sorted_indices[div(length(sorted_indices) + 1, 2)]
    
    (
        year_min = group.year[min_idx],
        year_med = group.year[median_idx],
        year_max = group.year[max_idx],
        min_extent = group.total_extents[min_idx],
        median_extent = group.total_extents[median_idx],
        max_extent = group.total_extents[max_idx]
    )
end

flood_years = [2011] #vcat(events.year_min,events.year_med, events.year_max)
one_shock = true
repeat_shocks = false

# Set up directories and logging
output_dir = joinpath(@__DIR__,"shap_$(ENV["SLURM_JOB_ID"])")
mkpath(output_dir)

# Set up logging files
run_log = joinpath(output_dir, "run_log.txt")
error_log = joinpath(output_dir, "error_log.txt")
progress_file = joinpath(output_dir, "progress.txt")

# Helper functions for logging
function log_info(msg)
    timestamp = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
    println(stdout, "INFO [$timestamp]: $msg")
    open(run_log, "a") do io
        println(io, "INFO [$timestamp]: $msg")
    end
    flush(stdout)
end

function log_error(msg)
    timestamp = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
    println(stderr, "ERROR [$timestamp]: $msg")
    open(error_log, "a") do io
        println(io, "ERROR [$timestamp]: $msg")
    end
    flush(stderr)
end

# Log start of process
log_info("Starting script with $(nworkers()) workers")
log_info("SLURM Job ID: $(get(ENV, "SLURM_JOB_ID", "unknown"))")
log_info("Results will be saved to: $output_dir")


# Set up additional parameters
add_params = OrderedDict(
    :pop_no=>[0,1,2,3,4,5,6,7,8,9,10],
    :seed=>collect(range(1000,1249))
)

# Extract parameter combinations and names
p_combs = collect((Tuple(row) for row in eachrow(calib_combs)))
output_params = collect(Symbol.(names(calib_combs)))
append!(output_params,collect(keys(add_params)))

# Generate all combinations outside of flood shock characteristics
mod_combs = [(c..., p, s) for (c, p, s) in Iterators.product(p_combs,values(add_params)...)];
param_matrix = stack([collect(tuple) for tuple in vec(mod_combs)])'
shap_param_df = DataFrame(param_matrix, output_params)

CSV.write(joinpath(output_dir,"param_runs_shap.csv"), shap_param_df)

append!(output_params, [:flood_event_year, :flood_repeat])

# Determine total combinations and chunk size
chunk_size = 39875  # Adjust based on memory requirements



for flood_shock in flood_years
    log_info("Starting flood shock year $flood_shock")
    
    #Create combinations with flood shock characteristics included
    combs = [(c..., p, s, flood_shock, repeat_shocks) for (c, p, s) in Iterators.product(p_combs,values(add_params)...)];
    n_combs = length(combs)
    n_chunks = ceil(Int, n_combs / chunk_size)

    # Write initial progress state
    open(progress_file, "w") do io
        println(io, "Starting flood shock year $flood_shock")
        println(io, "$(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))")
    end

    # Set up data file 
    filename = joinpath(output_dir,"$(flood_shock)_abm_data_shap_$(ENV["SLURM_JOB_ID"].h5")
    n_years = 40
    n_agents = 755

    h5open(filename, "w") do file
        # Create datasets with chunking for efficient I/O
        chunk_size = (1, n_agents, n_years, 1)
            
        # Main data array: (runs, agents, years, variables)
        create_dataset(file, "pop_data", Float32, (n_combs, n_agents, n_years, 3),
                        chunk=chunk_size, deflate=9, shuffle=true
        )

        create_dataset(file, "price_data", Float32, (n_combs, n_agents, n_years, 6),
                        chunk=chunk_size, deflate=9, shuffle=true
        )  
        # Metadata
        write(file, "pop_vars", string.(shap_adata[3:5]))
        write(file, "historical flood year", flood_shock)
        write(file, "n_runs", n_combs)
        write(file, "n_agents", n_agents)
        write(file, "GEOID", unique(phil_bg.GEOID)) 
        write(file, "n_years", n_years)
        write(file, "price_vars", string.(shap_adata[6:end]))    
    end

    # Set up pop shares during shock data file 
    pop_share_file = joinpath(output_dir,"$(flood_shock)_pop_share_data_shap_$(ENV["SLURM_JOB_ID"].h5")

    h5open(pop_share_file, "w") do file
        # Create datasets with chunking for efficient I/O
        chunk_size = (1,9,4)
            
        # Main data array: (runs, agents, years, variables)
        create_dataset(file, "pop_share_data", Float32, (n_combs, n_agents*9, 4),
                        chunk=chunk_size, deflate=9, shuffle=true
        )
    
        # Metadata
        write(file, "categories", ["low", "middle","high"])
        write(file, "column_names", ["GEOID", "agent group","housing cat", "count"])
        write(file, "n_runs", n_combs)
        write(file, "n_agents", n_agents)
        write(file, "GEOID", unique(phil_bg.GEOID)) 
        write(file, "n_years", n_years)
        write(file, "n_cat_combo", 9)
            
    end

   
    

    log_info("Processing $n_combs parameter combinations in $n_chunks chunks")

    for chunk_idx in 1:n_chunks
        # Get subset of combinations for this chunk
        start_idx = (chunk_idx - 1) * chunk_size + 1
        end_idx = min(chunk_idx * chunk_size, n_combs)
        chunk_combs = combs[start_idx:end_idx]
        
        log_info("Starting chunk $chunk_idx of $n_chunks (combinations $start_idx to $end_idx)")

        # Create a channel to stream results
        result_channel = RemoteChannel(() -> Channel{Tuple{Int, DataFrame, DataFrame}}(nworkers() * 2))

        # Async task to save results as they arrive
        save_task = @async begin
            valid_count = 0
            invalid_count = 0
            
            while true
                result = take!(result_channel)
                if isnothing(result)  # Sentinel value to stop
                    break
                end
                
                try
                    idx, sim_df, pop_df = result
                    save_model_data!(filename, idx, sim_df, n_agents, n_years)
                    save_pop_share_data!(pop_share_file, idx, pop_df)
                    valid_count += 1
                catch e
                    log_error("Error saving result $idx: $(sprint(showerror, e))")
                    invalid_count += 1
                end
            end
            
            (valid_count, invalid_count)
        end

        # Run simulations and stream results
        @sync begin
            @async begin
                # Set up progress meter
                progress = ProgressMeter.Progress(
                    length(chunk_combs); 
                    desc="Chunk $chunk_idx/$n_chunks: ",
                    enabled=true,
                    output=stderr,  # ProgressMeter outputs to stderr by default
                    dt=5.0  # Update every 5 seconds
                )
                
                ProgressMeter.progress_pmap(enumerate(chunk_combs); progress) do (i, comb)
                    try
                        result = run_single(comb, output_params, PhilPopABM; adata=shap_adata, n=39, shock=one_shock)
                        put!(result_channel, (start_idx + i - 1, result[1], result[2]))
                    catch e
                        worker_id = myid()
                        error_message = sprint(showerror, e, catch_backtrace())
                        println(stderr, "ERROR [Worker $worker_id]: $error_message")
                        put!(result_channel, (start_idx + i - 1, DataFrame(), DataFrame()))
                    end
                end
                put!(result_channel, nothing)  # Signal completion
            end
        end

        valid_count, invalid_count = fetch(save_task)
        log_info("Chunk $chunk_idx processed: $valid_count valid results, $invalid_count invalid results")
        
        # Close the channel to free resources
        close(result_channel)
        # Force garbage collection after processing chunk
        GC.gc()

        # Update progress file
        open(progress_file, "a") do io
            println(io, "Finished processing chunk $chunk_idx at $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))")
        end

        log_info("Flood shock year $flood_shock processed")
        open(joinpath(output_dir, "status.txt"), "w") do io
            println(io, "Flood shock year $flood_shock processed")
        end

        # Update status file for monitoring
        open(joinpath(output_dir, "status.txt"), "w") do io
            println(io, "Flood shock year $flood_shock processed")
            println(io, "Completed $chunk_idx of $n_chunks chunks")
            println(io, "Last update: $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))")
            println(io, "Valid results in last chunk: $valid_count")
            println(io, "Invalid results in last chunk: $invalid_count")
        end
    end
    #= Process each chunk
    for chunk_idx in 1:n_chunks
        # Get subset of combinations for this chunk
        start_idx = (chunk_idx - 1) * chunk_size + 1
        end_idx = min(chunk_idx * chunk_size, n_combs)
        chunk_combs = combs[start_idx:end_idx]
        
        log_info("Starting chunk $chunk_idx of $n_chunks (combinations $start_idx to $end_idx)")
        
        # Write initial progress state
        open(progress_file, "w") do io
            println(io, "Starting chunk $chunk_idx of $n_chunks")
            println(io, "$(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))")
        end
        
        # Set up progress meter
        progress = ProgressMeter.Progress(
            length(chunk_combs); 
            desc="Chunk $chunk_idx/$n_chunks: ",
            enabled=true,
            output=stderr,  # ProgressMeter outputs to stderr by default
            dt=5.0  # Update every 5 seconds
        )
        
        # Run simulations for this chunk
        chunk_results = try
            ProgressMeter.progress_pmap(chunk_combs; progress) do comb
                try
                    # Run Simulation
                    run_single(comb, output_params, PhilPopABM; adata=shap_adata, n=39, shock=one_shock)
                catch e
                    # Log worker errors but don't fail the whole chunk
                    worker_id = myid()
                    error_message = sprint(showerror, e, catch_backtrace())
                    # Cannot directly call log_error from workers, so we print to stderr
                    println(stderr, "ERROR [Worker $worker_id]: $error_message")
                    return (DataFrame(), DataFrame())  # Return empty dataframe on error
                end
            end
        catch e
            # Log main process errors
            log_error("Error in main process for chunk $chunk_idx: $(sprint(showerror, e, catch_backtrace()))")
            Tuple{DataFrame, DataFrame}[]
        end
        =#       
        
        
        # Check if we got any results
        #if isempty(chunk_results)
        #    log_error("No valid results for chunk $chunk_idx, skipping")
        #    continue
        #end
        #=
        # Process results immediately
        log_info("Processing results for chunk $chunk_idx")
        
        # Count valid and invalid results
        valid_count = 0
        invalid_count = 0
        
        for (i, result) in enumerate(chunk_results)
            try
                sim_df, pop_df = result
                save_model_data!(filename, (start_idx+i)-1, sim_df, n_agents, n_years)
                save_pop_share_data!(pop_share_file, (start_idx+i)-1, pop_df)
                valid_count += 1
            catch e
                log_error("Error processing result $i in chunk $chunk_idx: $(sprint(showerror, e))")
                invalid_count += 1
            end
        end
        =#
        
        
        
    
    #end

    
    
end


log_info("Script completed. All results saved to $output_dir")
