#!/bin/bash
#SBATCH -J comp
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=8G
#SBATCH -t 12:00:00
#SBATCH -o logs/comp_%j.log
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
module purge
module use "$HOME/modulefiles"
module load freefem/4.15 sqlite/3.45.0
mkdir -p logs
bash ./run_comp.sh ${AMP:-10} ${TEND:-10} ${TSCALE:-1} ${TAG:-auto} ${SHEAR:-1}
