#!/bin/bash
# time_convergence_dt8.sh -- one more point on the time-step ladder: dt/8 (0.00125) next to the existing dt, dt/2, dt/4 runs of time_convergence.sh.
# The relaxation checkpoint already exists (shared, made at dt 0.01), so nothing is relaxed again. Submits a chain of run.sbatch links and post.sbatch.
# usage (repo root, login node):  bash time_convergence_dt8.sh        env: TEND (5: dt/8 costs 8x the production dt, so it only runs to 5 ns; the common frames
#   of the ladder stop there), NBASE (5 links of 24 h at DT0, x8 here; links after the end exit at once), RUNTAG (kt3pcgtc), NP/MEM, HWCONSTRAINT/HWEXCLUDE, SKIPVERIFY=1
set -euo pipefail
cd "$(dirname "$0")"
export MODE=mpi ANGLE=90 TEND=${TEND:-5} RUNTAG=${RUNTAG:-kt3pcgtc}
export HOLD=${HOLD:-2.0} RAMP=${RAMP:-0.5} VTKEVERY=${VTKEVERY:-10}
export NP=${NP:-16} MEM=${MEM:-48G}
DT0=${DT0:-0.01}; NBASE=${NBASE:-5}
RES="-n $NP --mem=$MEM"
[ -n "${HWCONSTRAINT:-}" ] && RES="$RES --constraint=$HWCONSTRAINT"
[ -n "${HWEXCLUDE:-}" ] && RES="$RES --exclude=$HWEXCLUDE"
mkdir -p logs
if [ "${SKIPVERIFY:-0}" != 1 ] && [ -f checks/verify_solver.sh ]; then bash checks/verify_solver.sh || { echo "FAILED: verify_solver.sh (SKIPVERIFY=1 to bypass if you know why)"; exit 1; }; fi
bash plugin/build_plugin.sh || exit 1
export DT=$(awk -v d=$DT0 'BEGIN{printf "%g", d/8}')
TAG=$(UPMS=400; source ./tag.sh >/dev/null && echo "$TAG")
N=$((NBASE * 8))
prev=$(sbatch --parsable --export=ALL -J tc3 $RES --array=5 run.sbatch); first=$prev
for k in $(seq 2 "$N"); do prev=$(sbatch --parsable --export=ALL --array=5 --dependency=afterany:$prev -J tc3 $RES run.sbatch); done
post=$(sbatch --parsable --export=ALL --array=5 --dependency=afterany:$prev post.sbatch)
echo "dt $DT: $TAG  runs $first .. $prev ($N links)  post $post"
echo "progress: tail -n 1 out_$TAG/*_log.csv | cut -d, -f1-3"
echo "when done, the three finest runs form a second ladder:  RUNS=\"out_a90_up400ms_dt0.005_v2_kt3pcgtc out_a90_up400ms_dt0.0025_v2_kt3pcgtc out_$TAG\" bash time_convergence_norms.sh"
