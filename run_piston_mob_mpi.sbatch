#!/bin/bash
#SBATCH -J pistonMPI
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem=24G
#SBATCH -t 24:00:00
#SBATCH --array=0-5
#SBATCH -o logs/pistonMPI_%A_%a.log
# MPI/PETSc version (cosserat2d_mob_mpi.edp, -n ranks on one node; tag suffix _mpi, see run_piston_mob_mpi.sh) of the mobile-density drag variant of the confined-piston ladder: variant A only, cosserat2d_mob.edp.
#   defaults: DRAGRHO=3 (evolved mobile density), VC=2 (relativistic drag, km/s), GND=curl, RECOV=2 (Ask recovery, CD=100), KWCEDIFF=0, ETAPHIMIN=1e-2, ETASTAR=1 -> out_pistM3v2GRcV_gc_ne_pm1e-2_es1_a45_up<m/s>ms/
#   tasks 0-5: u_p = 25 50 100 150 250 400 m/s
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
echo "task $ID: variant A (dragrho=${DRAGRHO:-3}, vc=${VC:-2}, gnd=${GND:-curl}${GNDCURL:+ gndcurl=$GNDCURL}${GNDPHI2:+ gndphi2=$GNDPHI2}, recov=${RECOV:-2}), piston u_p = $UP m/s, tend ${TEND:-30} ns"
bash ./run_piston_mob_mpi.sh $UP ${TEND:-30} ${RAMP:-0.5} ${HOLD:-2.0} ${VTKEVERY:-10} A
