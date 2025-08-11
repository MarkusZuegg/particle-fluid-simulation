#!/bin/bash
#SBATCH --job-name=2D-sph-simulation
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=./logs/x%j.out
#SBATCH --error=./logs/x%j.err
#SBATCH --gpus-per-task=0

time ./sph $SLURM_JOB_ID
