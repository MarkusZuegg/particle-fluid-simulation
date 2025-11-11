#!/bin/bash
set -oue pipefail

mkdir -p ./logs

sbatch <<EOM
#!/bin/bash
#SBATCH --job-name=SPH
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=./logs/x%j.out
#SBATCH --error=./logs/x%j.err
#SBATCH --gpus-per-task=0
#SBATCH --time=00:10:00

time ./sph \$SLURM_JOB_ID
EOM
