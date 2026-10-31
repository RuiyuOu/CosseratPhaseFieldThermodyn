#!/bin/bash
# Full variant-A piston ladder (u_p = 25 50 100 150 250 400 m/s) with the current cosserat2d_mob.edp defaults,
# half time step and 10-step tangent refresh, followed automatically by the post-processing array.
#   physics defaults: DRAGRHO=3 VC=2 RECOV=2 CD=100 EOS=vinet VOLTAN=1 DISSIP=1 CDSTORE=1 RELAXFIX=1 ONFAIL=retry CKROLL=20
#   outputs: out_pistM3v2RcV_dt0.005_nr10_up<speed>ms/   animations: deliver_pistM3v2RcV_dt0.005_nr10_up<speed>ms/
# usage: bash run_all_A.sh [N_chained_arrays=20]   (full Newton + consistent tangent is ~10x slower per step than the old modified Newton)
#   each array job has 24 h; a level that is complete, hit TSTOP, stalled or stopped on a non-converged step exits at once,
#   so extra links cost nothing. Override any default on the command line, e.g.  EOS=svk bash run_all_A.sh
set -euo pipefail
cd "$(dirname "$0")"
N=${1:-20}
export DT=${DT:-0.005} NREF=${NREF:-10}
export TEND RAMP HOLD VTKEVERY TSTOP GNDPHI2 RECOV CA CD DIAG DRAGRHO VC RHOMSAT C2M RHOFCUT RHOMMIN RELAXFIX \
       SPINMAX SPINBRANCH ONFAIL MAXITRETRY CKROLL EXACTEXP DISSIP CDSTORE SOFTMIN EOS KP K0 VOLTAN NEWTON TOLRES RESFLOOR LSMAX MAXIT NLOCAL TANGENT POLISH FDK FDQ LINSOLVER LINTOL PRECEVERY TANEVERY SLIPPRED RUNTAG
[ -s out_A/A_dyn_relaxed_ck.txt ] || { echo "missing out_A/A_dyn_relaxed_ck.txt (run the variant-A relaxation first)"; exit 1; }
mkdir -p logs
prev=$(sbatch --parsable --export=ALL run_piston_mob.sbatch)
echo "array 1/$N: $prev"
for i in $(seq 2 "$N"); do
  prev=$(sbatch --parsable --export=ALL --dependency=afterany:$prev run_piston_mob.sbatch)
  echo "array $i/$N: $prev"
done
post=$(sbatch --parsable --export=ALL --dependency=afterany:$prev run_post_piston_mob.sbatch)
echo "post-processing array: $post (starts after $prev)"
echo "outputs: out_pistM3v2RcV_dt${DT}_nr${NREF}_up<speed>ms/   (check progress: tail -1 out_*/*_log.csv)"
