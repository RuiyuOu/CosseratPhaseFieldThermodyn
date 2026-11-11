#!/bin/bash
# Post-process every angle now (frames written so far), then the grain summary once all post jobs have finished.
# usage: bash post_all.sh
# env: ANGLES ("0 45 90"), ARRAY (0-5 = 25 50 100 150 250 400 m/s), FORCE=1 (redraw all animation frames),
#      MODE, RUNTAG, DT as for the runs, and the post variables of post.sbatch (TAUWIN0, TAUWIN1, TAUEVERY, GRAD, ...).
set -euo pipefail
cd "$(dirname "$0")"
export MODE=${MODE:-mpi} ANGLES="${ANGLES:-0 45 90}"
ARRAY=${ARRAY:-0-5}
mkdir -p logs
ids=""
for A in $ANGLES; do
  id=$(ANGLE=$A sbatch --parsable --export=ALL --array=$ARRAY post.sbatch)
  echo "angle $A: post-processing array $id"
  ids="$ids:$id"
done
echo "grain summary: $(sbatch --parsable --export=ALL --dependency=afterany$ids grain_summary.sbatch) (starts after$ids)"
UPMS=""; source ./tag.sh
echo "outputs: deliver_a<angle>_up<speed>ms${SUFFIX}/, grain_summary${SUFFIX}/"
