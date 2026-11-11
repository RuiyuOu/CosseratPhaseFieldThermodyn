#!/bin/bash
# One angle, all speeds: relaxation job, N chained 24 h arrays, post-processing.   usage: ANGLE=0 bash submit_speeds.sh [N]
# env: MODE (mpi: N 10, NP 16 ranks, MEM 48G | serial: N 20, 1 core, 8G, DT 0.005), ANGLE (45), ARRAY (0-5 = 25 50 100 150 250 400 m/s),
#      TEND RAMP HOLD VTKEVERY RUNTAG DT, physics (tag.sh), RELAXDEP (id of an already submitted relaxation job).
set -euo pipefail
cd "$(dirname "$0")"
export MODE=${MODE:-mpi} ANGLE=${ANGLE:-45}
ARRAY=${ARRAY:-0-5}
if [ "$MODE" = mpi ]; then
  N=${1:-10}; export NP=${NP:-16} MEM=${MEM:-48G}; RES="-n $NP --mem=$MEM"
  [ -s cosserat_kwc_mpi.edp ] || python3 mpi/make_mpi.py cosserat_kwc.edp cosserat_kwc_mpi.edp
else
  N=${1:-20}; export DT=${DT:-0.005}; RES="-n 1 --mem=${MEM:-8G}"
fi
mkdir -p logs
# relaxation stage: exits at once if relax_checkpoint.txt exists with the same parameters; the first array waits for it (afterok)
if [ -z "${RELAXDEP:-}" ]; then RELAXDEP=$(sbatch --parsable --export=ALL $RES run_relax.sbatch); echo "relaxation stage: $RELAXDEP"; fi
D0="afterok:$RELAXDEP"; [ -n "${SBATCH_DEPENDENCY:-}" ] && D0="$D0,$SBATCH_DEPENDENCY"
prev=$(sbatch --parsable --export=ALL --array=$ARRAY --dependency=$D0 $RES run.sbatch)
echo "array 1/$N: $prev"
for i in $(seq 2 "$N"); do
  prev=$(sbatch --parsable --export=ALL --array=$ARRAY --dependency=afterany:$prev $RES run.sbatch)
  echo "array $i/$N: $prev"
done
post=$(sbatch --parsable --export=ALL --array=$ARRAY --dependency=afterany:$prev post.sbatch)
echo "post-processing array: $post (starts after $prev)"
UPMS=""; source ./tag.sh
echo "outputs: out_a${ANGLE}_up<speed>ms${SUFFIX}/   (progress: tail -qn 1 out_a${ANGLE}_up*ms${SUFFIX}/*_log.csv)"
