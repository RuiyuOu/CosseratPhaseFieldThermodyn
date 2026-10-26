#!/bin/bash
#SBATCH -J piston
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=8G
#SBATCH -t 24:00:00
#SBATCH --array=0-5
#SBATCH -o logs/piston_%A_%a.log
# Confined-piston velocity ladder (cold-spray analogue), 30 ns each, all six levels in parallel (one array task each):
#   task: 0   1   2    3    4    5
#   u_p:  25  50  100  150  250  400 m/s   (impact v = 2 u_p; shock pressure ~ 0.9 1.8 3.7 5.6 9.7 16 GPa, Cu fit)
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
module purge
module use "$HOME/modulefiles"
module load freefem/4.15 sqlite/3.45.0
mkdir -p logs
UPS=(25 50 100 150 250 400)
UP=${UPS[${SLURM_ARRAY_TASK_ID:-0}]}
echo "task ${SLURM_ARRAY_TASK_ID:-0}: piston u_p = $UP m/s, tend ${TEND:-30} ns"
bash ./run_piston.sh $UP ${TEND:-30} ${RAMP:-0.5} ${HOLD:-2.0} ${VTKEVERY:-10}
