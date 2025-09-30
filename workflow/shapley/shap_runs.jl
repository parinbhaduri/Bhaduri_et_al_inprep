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
    using HDF5
end


@everywhere include(joinpath(dirname(@__DIR__),"src/data_collect.jl"))
@everywhere include("shap_functions.jl")




#Load calibrated parameter combinations
param_path = joinpath(dirname(@__DIR__),"calibration","data/param_comb_final_135842_mean_thresh_5_ens_250.csv")
calib_combs = DataFrame(CSV.File(param_path))[:,1:13]



# Set up directories and logging
output_dir = joinpath(@__DIR__,"data/shap_$(ENV["SLURM_JOB_ID"])")
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
# Generate all combinations
combs = [(c..., p, s) for (c, p, s) in Iterators.product(p_combs,values(add_params)...)]

# Determine total combinations and chunk size
chunk_size = 30250  # Adjust based on memory requirements
n_combs = length(combs)
n_chunks = ceil(Int, n_combs / chunk_size)

log_info("Processing $n_combs parameter combinations in $n_chunks chunks")

# Set up data file 
filename = joinpath(output_dir,"abm_data_$(ENV["SLURM_JOB_ID"]).h5") 
n_years = 40
var_names = string.(shap_adata[2:end])
n_vars = length(var_names)
n_agents = 755

h5open(filename, "w") do file
    # Create datasets with chunking for efficient I/O
    chunk_size = (min(100, n_combs), n_agents, min(10, n_years), n_vars)
        
    # Main data array: (runs, agents, years, variables)
    create_dataset(file, "data", Float64, (n_combs, n_agents, n_years, n_vars),
                      chunk=chunk_size, compress=3)
        
    # Metadata
    write(file, "variable_names", var_names)
    write(file, "n_runs", n_combs)
    write(file, "n_agents", n_agents)
    write(file, "GEOID", unique(phil_bg.GEOID)) 
    write(file, "n_years", n_years)
    write(file, "n_variables", n_vars)
        
end

# Initialize tracking variables
#adf_cols = nothing
#mdf_cols = nothing
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
                # Run Simulation
                run_single(comb, output_params, PhilPopABM; adata=shap_adata, n=39)
            catch e
                # Log worker errors but don't fail the whole chunk
                worker_id = myid()
                error_message = sprint(showerror, e, catch_backtrace())
                # Cannot directly call log_error from workers, so we print to stderr
                println(stderr, "ERROR [Worker $worker_id]: $error_message")
                return (DataFrame())  # Return empty dataframe on error
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
    
    for (i, result) in enumerate(chunk_results)
        try
            save_model_data!(filename, (start_idx+i)-1, result, n_agents, n_years, n_vars)
            valid_count += 1
        catch e
            log_error("Error processing result $i in chunk $chunk_idx: $(sprint(showerror, e))")
            invalid_count += 1
        end
    end
    
    log_info("Chunk $chunk_idx processed: $valid_count valid results, $invalid_count invalid results")
    
    # Clear memory
    chunk_results = nothing
    GC.gc()
    
    # Update status file for monitoring
    open(joinpath(output_dir, "status.txt"), "w") do io
        println(io, "Completed $chunk_idx of $n_chunks chunks")
        println(io, "Last update: $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))")
        println(io, "Valid results in last chunk: $valid_count")
        println(io, "Invalid results in last chunk: $invalid_count")
    end
   
end

log_info("Script completed. All results saved to $output_dir")