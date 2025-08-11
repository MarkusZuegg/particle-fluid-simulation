#!/bin/bash
#SBATCH --job-name=2D-sph-simulation
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=x%j.out
#SBATCH --error=x%j.err

srun ./sph $SLURM_JOB_ID
