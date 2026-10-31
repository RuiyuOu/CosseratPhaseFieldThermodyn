#!/bin/bash
#SBATCH -J pistonM
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=8G
#SBATCH -t 24:00:00
#SBATCH --array=0-5
#SBATCH -o logs/pistonM_%A_%a.log
# Mobile-density drag variant of the confined-piston ladder: variant A only, cosserat2d_mob.edp.
#   defaults: DRAGRHO=3 (evolved mobile density), VC=2 (relativistic drag, km/s), GNDPHI2=1, RECOV=2 (Ask recovery, CD=100) -> out_pistM3v2RcV_up<m/s>ms/
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
echo "task $ID: variant A (dragrho=${DRAGRHO:-3}, vc=${VC:-2}, gndphi2=${GNDPHI2:-1}, recov=${RECOV:-2}), piston u_p = $UP m/s, tend ${TEND:-30} ns"
bash ./run_piston_mob.sh $UP ${TEND:-30} ${RAMP:-0.5} ${HOLD:-2.0} ${VTKEVERY:-10} A
