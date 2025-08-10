#!/bin/bash
#SBATCH --job-name=2D-sph-simulation
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

time ./mile0 {JOBID}
