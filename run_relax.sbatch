#!/bin/bash
#SBATCH -J relax
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mem=48G
#SBATCH -t 24:00:00
#SBATCH -o logs/relax_%j.log
# Relaxation stage before a ladder (see run_relax.sh): exits at once if relax_checkpoint.txt already exists with the same parameters.
# The submit scripts pass -n and --mem; MODE=serial (exported) runs the serial solver instead.
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
module purge
module use "$HOME/modulefiles"
module load freefem/4.15 sqlite/3.45.0
# (load here the MPI module FreeFem++-mpi was built with, if the freefem module does not load it)
mkdir -p logs
bash ./run_relax.sh
