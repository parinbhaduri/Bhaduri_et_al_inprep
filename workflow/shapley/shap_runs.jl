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

@everywhere begin 
    include(joinpath(dirname(@__DIR__),"src","data_include.jl"))
    include(joinpath(dirname(@__DIR__),"src","functions.jl"))
    include(joinpath(dirname(@__DIR__),"src","data_collect.jl"))
end


#Load calibrated parameter combinations
param_path = joinpath(dirname(@__DIR__),"calibration","data/param_comb_final_mean_thresh_6_ens_250.csv")
calib_combs = DataFrame(CSV.File(param_path))[:,1:14]
#=
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
=#
flood_years = [1991,2018,1981] #vcat(events.year_min,events.year_med, events.year_max)
one_shock = true
repeat_shocks = false

# Set up directories and logging
output_dir = joinpath(@__DIR__,"data","shap_$(ENV["SLURM_JOB_ID"])")
mkpath(output_dir)

data_dir = joinpath(@__DIR__,"data","shap_runs")
mkpath(data_dir)

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

append!(output_params, [:flood_event_year, :flood_shock])

# Determine total combinations and chunk size
chunk_size = 39875  # Adjust based on memory requirements



for flood_shock in flood_years
    log_info("Starting flood shock year $flood_shock")
    
    #Create combinations with flood shock characteristics included
    combs = [(c..., p, s, flood_shock, one_shock) for (c, p, s) in Iterators.product(p_combs,values(add_params)...)];
    n_combs = length(combs)
    n_chunks = ceil(Int, n_combs / chunk_size)

    # Write initial progress state
    open(progress_file, "w") do io
        println(io, "Starting flood shock year $flood_shock")
        println(io, "$(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))")
    end

    # Set up data file 
    filename = joinpath(data_dir,"$(flood_shock)_abm_data.h5")
    n_years = 20
    n_agents = 755
    
    h5open(filename, "w") do file
        # Create datasets with chunking for efficient I/O
        chunk_size = (1, n_agents, n_years+1, 1)
            
        # Main data array: (runs, agents, years, variables)
        create_dataset(file, "pop_data", Float32, (n_combs, n_agents, n_years+1, 6),
                        chunk=chunk_size, deflate=9, shuffle=true
        )

        create_dataset(file, "price_data", Float32, (n_combs, n_agents, n_years+1, 6),
                        chunk=chunk_size, deflate=9, shuffle=true
        )  
        # Metadata
        write(file, "pop_vars", string.(shap_adata[3:8]))
        write(file, "historical flood year", flood_shock)
        write(file, "n_runs", n_combs)
        write(file, "n_agents", n_agents)
        write(file, "GEOID", unique(phil_bg.GEOID)) 
        write(file, "n_years", n_years)
        write(file, "price_vars", string.(shap_adata[9:end]))    
    end
    
    # Set up pop shares during shock data file 
    pop_share_file = joinpath(data_dir,"$(flood_shock)_pop_share_data.h5")

    h5open(pop_share_file, "w") do file
        # Create datasets with chunking for efficient I/O
        chunk_size = (1,9,7) 
            
        # Main data array: (runs, agents, years, variables)
        create_dataset(file, "pop_share_data", Float64, (n_combs, n_agents*9, 7), 
                        chunk=chunk_size, deflate=9, shuffle=true
        )
    
        # Metadata
        write(file, "categories", ["low", "middle","high"])
        write(file, "column_names", ["GEOID", "agent group","housing cat", "agent count", "income", "household count", "pop. count"])
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
        result_channel = RemoteChannel(() -> Channel{Tuple{Int, DataFrame, DataFrame}}(nworkers() * 2)) # 

        # Async task to save results as they arrive
        save_task = @async begin
            valid_count = 0
            invalid_count = 0
            processed_count = 0
    
            while processed_count < length(chunk_combs)
                try
                    idx, sim_df, pop_df = take!(result_channel)
                    #idx, sim_df, pop_df = results
                    try  
                        save_model_data!(filename, idx, sim_df, n_agents, n_years+1)
                        save_pop_share_data!(pop_share_file, idx, pop_df)
                        valid_count += 1
                    catch e
                        log_error("Error saving result $idx: $(sprint(showerror, e))")
                        invalid_count += 1
                    end

                    processed_count += 1

                catch e
                    log_error("Error taking from channel: $(sprint(showerror, e))")
                    break
                end
            end
            
            (valid_count, invalid_count)
        end

        # Run simulations and stream results
        
        # Set up progress meter
        progress = ProgressMeter.Progress(
            length(chunk_combs); 
            desc="Chunk $chunk_idx/$n_chunks: ",
            enabled=true,
            output=stderr,  # ProgressMeter outputs to stderr by default
            dt=5.0  # Update every 5 seconds
        )
        # Run simulations
        ProgressMeter.progress_pmap(enumerate(chunk_combs); progress) do (i, comb)
            try
                sim_df, pop_df = run_shap_single(comb, output_params, PhilPopABM; adata=shap_adata, n=20, shock=one_shock)
                put!(result_channel, (start_idx + i - 1, sim_df, pop_df))
            catch e
                worker_id = myid()
                error_message = sprint(showerror, e, catch_backtrace())
                println(stderr, "ERROR [Worker $worker_id]: $error_message")
                put!(result_channel, (start_idx + i - 1, DataFrame(), DataFrame()))
            end
        end

        valid_count, invalid_count = fetch(save_task)
        log_info("Chunk $chunk_idx processed: $valid_count valid results, $invalid_count invalid results")
        
        # Now safe to close the channel
        close(result_channel)
        # Force garbage collection after processing chunk
        GC.gc()

        # Update progress file
        open(progress_file, "a") do io
            println(io, "Finished processing chunk $chunk_idx at $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))")
        end

        log_info("Flood shock year $flood_shock processed (chunk $chunk_idx of $n_chunks)")

        # Update status file for monitoring
        open(joinpath(output_dir, "status.txt"), "w") do io
            println(io, "Flood shock year $flood_shock processed (chunk $chunk_idx of $n_chunks)")
            println(io, "Completed $chunk_idx of $n_chunks chunks")
            println(io, "Last update: $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))")
            println(io, "Valid results in last chunk: $valid_count")
            println(io, "Invalid results in last chunk: $invalid_count")
        end
    end
end


log_info("Script completed. All results saved to $output_dir")
