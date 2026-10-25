#!/bin/bash
# Submits the 7-level ladder as an array (all levels in parallel), then N-1 follow-up arrays (afterany) so levels that did not
# finish within 12 h resume from their checkpoints; finished levels exit at once; a level that STALLED (e.g. a NaN blow-up) is skipped.
# usage: bash submit_ladder.sh [N_arrays=4]     optional env: TEND (ns, default 10)
N=${1:-4}; mkdir -p logs; export TEND
prev=$(sbatch --parsable --export=ALL run_ladder.sbatch)
for i in $(seq 2 $N); do prev=$(sbatch --parsable --export=ALL --dependency=afterany:$prev run_ladder.sbatch); done
echo "ladder submitted; last array job $prev"
echo "outputs: out_lad_P<MPa>MPa_S<MPa>MPa/<same>_s*.vtu   (check each *_log.csv; a STALLED file marks a failed level)"
