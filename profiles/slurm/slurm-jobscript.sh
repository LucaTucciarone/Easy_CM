#!/bin/bash
# =============================================================================
# EasyCM -- SLURM job wrapper
# =============================================================================

#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=30G
#SBATCH --time=4:00:00
#SBATCH --account=csd854
#SBATCH --partition=condo
#SBATCH --qos=condo

# Activate conda environment
source /tscc/nfs/home/ltucciarone/miniconda3/etc/profile.d/conda.sh
conda activate EasyCM

# Run the Snakemake job
{exec_job}
