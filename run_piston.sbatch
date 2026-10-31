#!/bin/bash
#SBATCH -J piston
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=8G
#SBATCH -t 24:00:00
#SBATCH --array=0-11
#SBATCH -o logs/piston_%A_%a.log
# Confined-piston velocity ladder (cold-spray analogue), 30 ns each, all 12 jobs in parallel (one array task each):
#   tasks 0-5  = variant A, tasks 6-11 = variant F, u_p = 25 50 100 150 250 400 m/s (A and F use the same speeds in that order)   (impact v = 2 u_p; shock pressure ~ 0.9 1.8 3.7 5.6 9.7 16 GPa, Cu fit)
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
module purge
module use "$HOME/modulefiles"
module load freefem/4.15 sqlite/3.45.0
mkdir -p logs
UPS=(25 50 100 150 250 400)
ID=${SLURM_ARRAY_TASK_ID:-0}; VAR=A; if [ $ID -ge 6 ]; then VAR=F; ID=$((ID-6)); fi
UP=${UPS[$ID]}
echo "task ${SLURM_ARRAY_TASK_ID:-0}: variant $VAR, piston u_p = $UP m/s, tend ${TEND:-30} ns"
bash ./run_piston.sh $UP ${TEND:-30} ${RAMP:-0.5} ${HOLD:-2.0} ${VTKEVERY:-10} $VAR
