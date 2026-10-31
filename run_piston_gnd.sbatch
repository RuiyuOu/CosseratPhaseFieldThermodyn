#!/bin/bash
#SBATCH -J pistonG
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=8G
#SBATCH -t 24:00:00
#SBATCH --array=0-5
#SBATCH -o logs/pistonG_%A_%a.log
# GND/recovery variant of the confined-piston ladder: variant A only, cosserat2d_gnd.edp.
#   defaults: GNDPHI2=0 (rho_G without phi^2 weight), RECOV=2 (Ask 2019 recovery) -> out_pistGR_up<m/s>ms/; GNDPHI2=1 -> pistR, RECOV=0 -> pistG
#   tasks 0-5: u_p = 25 50 100 150 250 400 m/s
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
module purge
module use "$HOME/modulefiles"
module load freefem/4.15 sqlite/3.45.0
mkdir -p logs
UPS=(25 50 100 150 250 400)
ID=${SLURM_ARRAY_TASK_ID:-0}
UP=${UPS[$ID]}
echo "task $ID: variant A (gndphi2=${GNDPHI2:-0}, recov=${RECOV:-2}), piston u_p = $UP m/s, tend ${TEND:-30} ns"
bash ./run_piston_gnd.sh $UP ${TEND:-30} ${RAMP:-0.5} ${HOLD:-2.0} ${VTKEVERY:-10} A
