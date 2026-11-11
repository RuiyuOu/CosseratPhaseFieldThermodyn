#!/bin/bash
# Angles 0 45 90 x speeds 25..400 m/s: one relaxation job, the speed ladder per angle (submit_speeds.sh), then the grain summary.
# usage: bash submit_ladder.sh [N chained 24 h arrays]
# env: MODE (mpi | serial), ANGLES ("0 45 90"), ARRAY, NP, MEM, TEND, RUNTAG, SERIAL_ANGLES=1 (one angle at a time), post and physics variables.
#      Restart-safe: running it again resumes unfinished runs and skips finished ones.
set -euo pipefail
cd "$(dirname "$0")"
export MODE=${MODE:-mpi} ARRAY=${ARRAY:-0-5}
if [ "$MODE" = mpi ]; then
  export NP=${NP:-16} MEM=${MEM:-48G}; RES="-n $NP --mem=$MEM"; WHAT="$NP ranks, $MEM"
  [ -s cosserat_kwc_mpi.edp ] || python3 mpi/make_mpi.py cosserat_kwc.edp cosserat_kwc_mpi.edp
else
  export DT=${DT:-0.005}; RES="-n 1 --mem=${MEM:-8G}"; WHAT="serial"
fi
mkdir -p logs
export RELAXDEP=$(sbatch --parsable --export=ALL $RES run_relax.sbatch); echo "relaxation stage: $RELAXDEP (the ladder waits for it)"
dep=""; alldeps=""
for A in ${ANGLES:-0 45 90}; do
  echo "== angle $A: ARRAY=$ARRAY, $WHAT"
  if [ "${SERIAL_ANGLES:-0}" = 1 ] && [ -n "$dep" ]; then export SBATCH_DEPENDENCY="afterany:$dep"; fi
  out=$(ANGLE=$A bash submit_speeds.sh "$@")
  echo "$out"
  dep=$(echo "$out" | sed -n 's/^post-processing array: \([0-9]*\).*/\1/p')
  [ -n "$dep" ] && alldeps="$alldeps:$dep"
done
export ANGLES="${ANGLES:-0 45 90}"
[ -n "$alldeps" ] && echo "grain summary: $(sbatch --parsable --export=ALL --dependency=afterany$alldeps grain_summary.sbatch)"
UPMS=""; source ./tag.sh
echo "progress:  tail -qn 1 out_a*_up*ms${SUFFIX}/*_log.csv     queue: squeue -u \$USER"
