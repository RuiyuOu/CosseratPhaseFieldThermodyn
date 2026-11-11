#!/bin/bash
#SBATCH -J impact
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mem=48G
#SBATCH -t 24:00:00
#SBATCH --array=0-5
#SBATCH -o logs/impact_%A_%a.log
# One speed per array task (0-5 = 25 50 100 150 250 400 m/s) -> run.sh.  submit_speeds.sh sets -n/--mem (serial: -n 1 --mem=8G).
# env: MODE (mpi | serial), ANGLE (45), TEND (30), RAMP (0.5), HOLD (2.0), VTKEVERY (10), RUNTAG, DT, physics variables of tag.sh
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
module purge
module use "$HOME/modulefiles"
module load freefem/4.15 sqlite/3.45.0
# (load here the MPI module FreeFem++-mpi was built with, if the freefem module does not load it)
mkdir -p logs
UPS=(25 50 100 150 250 400)
ID=${SLURM_ARRAY_TASK_ID:-0}
UP=${UPS[$ID]}
echo "task $ID: plate $UP m/s at ${ANGLE:-45} deg, tend ${TEND:-30} ns, ${MODE:-mpi}"
bash ./run.sh $UP ${TEND:-30} ${RAMP:-0.5} ${HOLD:-2.0} ${VTKEVERY:-10}
