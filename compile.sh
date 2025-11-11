#!/bin/bash
set -euo pipefail

module load cuda/11.1

#I would have expected the module loads to add these, but apparently not
export PATH=/opt/local/stow/cuda-11.1/bin:$PATH

make clean
make all