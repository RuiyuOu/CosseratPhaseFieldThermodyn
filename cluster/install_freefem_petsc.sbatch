#!/bin/bash
#SBATCH -J ff-install
#SBATCH -p small
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 32
#SBATCH --mem=64G
#SBATCH -t 8:00:00
#SBATCH -o logs/ff_install_%j.log
# Heavy stages of the FreeFEM + PETSc build on a compute node (run fetch, configure and getall on the login node first).
# usage (repo root):  mkdir -p logs; sbatch cluster/install_freefem_petsc.sbatch          env: STAGES ("petsc build install module"), as install_freefem_petsc.sh
set -uo pipefail
cd "$SLURM_SUBMIT_DIR"
bash cluster/install_freefem_petsc.sh ${STAGES:-petsc build install module}
