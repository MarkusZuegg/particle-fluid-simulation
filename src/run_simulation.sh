#!/bin/bash
set -oue pipefail

if [[ -z "$1" ]]; then
    echo "Usage: $0 <cpu|gpu>"
    exit 1
fi

TARGET_EXEC="sph_$1"
REQUEST_GPUS=0

if [[ "$1" == "gpu" ]]; then
    REQUEST_GPUS=1
fi

mkdir -p ./logs

sbatch <<EOM
#!/bin/bash
#SBATCH --job-name=${TARGET_EXEC}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=./logs/%x%j.out
#SBATCH --error=./logs/%x%j.err
#SBATCH --gpus-per-task=${REQUEST_GPUS}
#SBATCH --time=00:10:00

module load cuda/11.1
export PATH=/opt/local/stow/cuda-11.1/bin:$PATH

time ./${TARGET_EXEC} \$SLURM_JOB_ID
EOM
