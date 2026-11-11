#!/bin/bash
#SBATCH -J adcompare
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mem=48G
#SBATCH -t 12:00:00
#SBATCH -o logs/adcompare_%j.log
# Stable (fd) vs AD tangent on the production set-up: 16 MPI ranks, same relaxed checkpoint, same plate (UPMS, ANGLE), TEND ns.
#   1. builds the AD variant, 2. -adcheck (serial: AD vs central differences on every element), 3. run.sh with TANGENT=fd and TANGENT=ad
#   (RUNTAG fdcmp / adcmp), one after the other on the same cores, 4. ad/compare_runs.py: per-step differences, iterations, wall time.
# usage (repo root, relax_checkpoint.txt present):  sbatch ad/ad_compare.sbatch      env: UPMS (400), ANGLE (45), TEND (1.0)
set -uo pipefail
cd "$SLURM_SUBMIT_DIR"
module purge
module use "$HOME/modulefiles"
module load freefem/4.15 sqlite/3.45.0
mkdir -p logs
UPMS=${UPMS:-400}; export ANGLE=${ANGLE:-45}; TEND=${TEND:-1.0}
bash ad/build_ad.sh || exit 1
echo "== adcheck (serial)"
mkdir -p ad/check && FreeFem++ -nw -v 0 cosserat_kwc_ad.edp -adcheck 1 -nrelax 0 -out ad/check/adcheck 2>&1 | grep "adcheck"
for T in fd ad; do
  echo "== TANGENT=$T"
  MODE=mpi TANGENT=$T RUNTAG=${T}cmp bash ./run.sh $UPMS $TEND 0.5 2.0 10
done
python3 ad/compare_runs.py out_a${ANGLE}_up${UPMS}ms_fdcmp out_a${ANGLE}_up${UPMS}ms_adcmp
