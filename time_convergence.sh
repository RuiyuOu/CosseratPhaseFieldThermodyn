#!/bin/bash
# time_convergence.sh -- time-step convergence test of the production set-up: 90 deg, 400 m/s, t_end 15 ns, at the production dt,
# dt/2 and dt/4.  Submits: one relaxation (shared: it always runs at dt 0.01), for each dt a chain of 24 h run.sbatch links
# (array task 5 = 400 m/s) and post.sbatch, then post/time_convergence.py once all three post jobs are done.
# usage (repo root, login node):  bash time_convergence.sh
# env:  RUNTAG (kt3pcgtc: its own run folders, so the production runs are untouched), DT0 (0.01, the production dt), NBASE (5
#       links of 24 h at DT0; dt/2 gets 2x, dt/4 4x; links after the end exit at once), TEND (15), HOLD/RAMP/VTKEVERY (production
#       defaults 2.0/0.5/10; frames keep their ns spacing at every dt), NP/MEM, HWCONSTRAINT/HWEXCLUDE/HWEXCL as submit_speeds.sh,
#       EVERY (compare fields every 0.5 ns), GRAINT ("3 6 9 12 15").
# output: out_a90_up400ms[_dt<dt>]_v2_<RUNTAG>/, deliver_..., and post/time_convergence/ (see post/time_convergence.py)
set -euo pipefail
cd "$(dirname "$0")"
export MODE=mpi ANGLE=90 TEND=${TEND:-15} RUNTAG=${RUNTAG:-kt3pcgtc}
export HOLD=${HOLD:-2.0} RAMP=${RAMP:-0.5} VTKEVERY=${VTKEVERY:-10}
export NP=${NP:-16} MEM=${MEM:-48G}
DT0=${DT0:-0.01}; NBASE=${NBASE:-5}
RES="-n $NP --mem=$MEM"
[ -n "${HWCONSTRAINT:-}" ] && RES="$RES --constraint=$HWCONSTRAINT"
[ -n "${HWEXCLUDE:-}" ] && RES="$RES --exclude=$HWEXCLUDE"
[ "${HWEXCL:-0}" = 1 ] && RES="$RES --exclusive"
mkdir -p logs post/time_convergence
if [ -x checks/verify_solver.sh ] || [ -f checks/verify_solver.sh ]; then
  bash checks/verify_solver.sh || { echo "FAILED: solver files differ from the frozen production set (checks/verify_solver.sh)"; exit 1; }
else
  echo "WARNING: checks/verify_solver.sh not found: cannot confirm the frozen solver"
fi
bash plugin/build_plugin.sh || exit 1
DTS=("$DT0" "$(awk -v d=$DT0 'BEGIN{printf "%g", d/2}')" "$(awk -v d=$DT0 'BEGIN{printf "%g", d/4}')")
RELAX=$(DT=$DT0 sbatch --parsable --export=ALL $RES run_relax.sbatch)
echo "relaxation: $RELAX"
POSTS=""; DIRS=""
for i in 0 1 2; do
  export DT=${DTS[$i]}
  N=$((NBASE * (1 << i)))
  TAG=$(UPMS=400; source ./tag.sh >/dev/null && echo "$TAG")
  prev=$(sbatch --parsable --export=ALL --array=5 --dependency=afterok:$RELAX -J tc$i $RES run.sbatch)
  first=$prev
  for k in $(seq 2 "$N"); do prev=$(sbatch --parsable --export=ALL --array=5 --dependency=afterany:$prev -J tc$i $RES run.sbatch); done
  post=$(sbatch --parsable --export=ALL --array=5 --dependency=afterany:$prev post.sbatch)
  echo "dt $DT: $TAG  runs $first .. $prev ($N links)  post $post"
  POSTS="$POSTS:$post"; DIRS="$DIRS out_$TAG"
done
unset DT
AN=$(sbatch --parsable -J tconv -p general -n 4 --mem=24G -t 08:00:00 -o logs/time_convergence_%j.log --dependency=afterany$POSTS \
  --wrap "source ~/dynpost_venv/bin/activate && NP=4 SYM=180 python3 post/time_convergence.py --runs$DIRS --out post/time_convergence --every ${EVERY:-0.5} --grain-times ${GRAINT:-3 6 9 12 15}")
echo "analysis: $AN (after$POSTS) -> post/time_convergence/  (log: logs/time_convergence_$AN.log)"
echo "progress: tail -qn 1$(for d in $DIRS; do printf ' %s/*_log.csv' "$d"; done) | cut -d, -f1-3"
