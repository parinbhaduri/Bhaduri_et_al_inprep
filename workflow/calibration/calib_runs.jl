#activate project environment
import Pkg
Pkg.activate(".")
Pkg.instantiate()


using Distributed, SlurmClusterManager


addprocs(SlurmManager())


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
end


@everywhere include(joinpath(dirname(@__DIR__),"src/data_collect.jl"))
@everywhere include("calib_functions.jl")




calib_params = OrderedDict(
    :seed=>collect(range(1000,1499)),
    :risk_averse=>[0.3,0.7], #0.5,
    :build_inc_perc=>[0.01, 0.02], #0.25,
    :price_inc_perc=>[0.01, 0.02], #
    :rhea_coef=>[0.65, 0.75],
    :base_move=>[0.01,0.03], #0.02,
    :prop_l=>[0.25, 0.75],
    :env_amen_l=>[0.25, 0.75],
    :prop_m=>[0.25, 0.75],
    :env_amen_m=>[0.25, 0.75],
    :prop_h=>[0.25, 0.75],
    :env_amen_h=>[0.25, 0.75],
    :penalty=>[0.25, 0.75],
    :flood_coefficient=>[0.25, 0.75]
)

chunk_size = 256000  # Adjust based on memory requirements

# Set up directories and logging
output_dir = joinpath(@__DIR__,"data/results_$(ENV["SLURM_JOB_ID"])")
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


# Extract parameter combinations and names
combs = collect(Iterators.product(values(calib_params)...))
output_params = collect(keys(calib_params))

# Determine total combinations and chunk size
n_combs = length(combs)
n_chunks = ceil(Int, n_combs / chunk_size)

log_info("Processing $n_combs parameter combinations in $n_chunks chunks")

# Initialize tracking variables
adf_cols = nothing
mdf_cols = nothing
first_run = true

# Process each chunk
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
                # Call your simulation function - REPLACE WITH YOUR ACTUAL FUNCTION CALL
                run_single(comb, output_params, PhilABM; adata=calib_adata, mdata=calib_mdata, n=39)
            catch e
                # Log worker errors but don't fail the whole chunk
                worker_id = myid()
                error_message = sprint(showerror, e, catch_backtrace())
                # Cannot directly call log_error from workers, so we print to stderr
                println(stderr, "ERROR [Worker $worker_id]: $error_message")
                return (DataFrame(), DataFrame())  # Return empty dataframes on error
            end
        end
    catch e
        # Log main process errors
        log_error("Error in main process for chunk $chunk_idx: $(sprint(showerror, e, catch_backtrace()))")
        Tuple{DataFrame, DataFrame}[]
    end
    
    # Update progress file
    open(progress_file, "a") do io
        println(io, "Finished processing chunk $chunk_idx at $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))")
    end
    
    # Check if we got any results
    if isempty(chunk_results)
        log_error("No valid results for chunk $chunk_idx, skipping")
        continue
    end
    
    # Process results immediately
    log_info("Processing results for chunk $chunk_idx")
    
    # Count valid and invalid results
    valid_count = 0
    invalid_count = 0
    
    adf_chunk = DataFrame()
    mdf_chunk = DataFrame()
    
    for (i, result) in enumerate(chunk_results)
        try
            df1, df2 = result
            if nrow(df1) > 0 && nrow(df2) > 0
                append!(adf_chunk, df1)
                append!(mdf_chunk, df2)
                valid_count += 1
            else
                invalid_count += 1
            end
        catch e
            log_error("Error processing result $i in chunk $chunk_idx: $(sprint(showerror, e))")
            invalid_count += 1
        end
    end
    
    log_info("Chunk $chunk_idx processed: $valid_count valid results, $invalid_count invalid results")
    
    # If first run with valid data, save column names
    if first_run && !isempty(adf_chunk) && !isempty(mdf_chunk)
        global adf_cols = names(adf_chunk)
        global mdf_cols = names(mdf_chunk)
        global first_run = false
        
        # Log column structure
        log_info("Agent dataframe columns: $(join(adf_cols, ", "))")
        log_info("Model dataframe columns: $(join(mdf_cols, ", "))")
    end
    
    # Save chunk results to CSV if we have valid data
    if !isempty(adf_chunk) && !isempty(mdf_chunk)
        chunk_filename_a = joinpath(output_dir, "agents_chunk_$(chunk_idx).csv")
        chunk_filename_m = joinpath(output_dir, "model_chunk_$(chunk_idx).csv")
        
        try
            CSV.write(chunk_filename_a, adf_chunk)
            CSV.write(chunk_filename_m, mdf_chunk)
            log_info("Chunk $chunk_idx saved successfully")
        catch e
            log_error("Error saving chunk $chunk_idx: $(sprint(showerror, e))")
        end
    else
        log_error("Chunk $chunk_idx had no valid data to save")
    end
    
    # Clear memory
    chunk_results = nothing
    adf_chunk = nothing
    mdf_chunk = nothing
    GC.gc()
    
    # Update status file for monitoring
    open(joinpath(output_dir, "status.txt"), "w") do io
        println(io, "Completed $chunk_idx of $n_chunks chunks")
        println(io, "Last update: $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))")
        println(io, "Valid results in last chunk: $valid_count")
        println(io, "Invalid results in last chunk: $invalid_count")
    end
end

log_info("All chunks processed. Creating final summary...")

# Create summary of available results
try
    chunk_files_a = filter(f -> startswith(f, "agents_chunk_") && endswith(f, ".csv"), readdir(output_dir))
    chunk_files_m = filter(f -> startswith(f, "model_chunk_") && endswith(f, ".csv"), readdir(output_dir))
    
    # Get row counts from the first row of each file to estimate total
    agent_row_counts = [
        try
            f = CSV.File(joinpath(output_dir, file), limit=1)
            length(Tables.rows(f))
        catch
            0
        end
        for file in chunk_files_a
    ]
    
    model_row_counts = [
        try
            f = CSV.File(joinpath(output_dir, file), limit=1)
            length(Tables.rows(f))
        catch
            0
        end
        for file in chunk_files_m
    ]
    
    open(joinpath(output_dir, "summary.txt"), "w") do io
        println(io, "Run completed at: $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))")
        println(io, "Total chunks: $n_chunks")
        println(io, "Agent data chunks saved: $(length(chunk_files_a))")
        println(io, "Model data chunks saved: $(length(chunk_files_m))")
        println(io, "Successfully processed chunks: $(length(filter(c -> c > 0, agent_row_counts)))")
        println(io, "Failed chunks: $(length(filter(c -> c == 0, agent_row_counts)))")
    end
    
    log_info("Summary created successfully")
catch e
    log_error("Error creating summary: $(sprint(showerror, e))")
end

log_info("Script completed. All results saved to $output_dir")