#!/bin/bash
#SBATCH --job-name=shap_outcome_drivers
#SBATCH --mail-user=pbb62@cornell.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=5
#SBATCH --mem-per-cpu=15G        # Increase memory per worker
#SBATCH --output=workflow/shapley/output_text.txt
#SBATCH --error=workflow/shapley/error_text.txt
#SBATCH --exclusive

# Print properties of job as submitted
echo "SLURM_JOB_ID = $SLURM_JOB_ID"
echo "SLURM_NTASKS = $SLURM_NTASKS"
echo "SLURM_NTASKS_PER_NODE = $SLURM_NTASKS_PER_NODE"
echo "SLURM_CPUS_PER_TASK = $SLURM_CPUS_PER_TASK"
echo "SLURM_JOB_NUM_NODES = $SLURM_JOB_NUM_NODES"

# Print properties of job as scheduled by Slurm
echo "SLURM_JOB_NODELIST = $SLURM_JOB_NODELIST"
echo "SLURM_TASKS_PER_NODE = $SLURM_TASKS_PER_NODE"
echo "SLURM_JOB_CPUS_PER_NODE = $SLURM_JOB_CPUS_PER_NODE"
echo "SLURM_CPUS_ON_NODE = $SLURM_CPUS_ON_NODE"


# Run the Julia code
julia +1.10 workflow/shapley/outcome_drivers.jl