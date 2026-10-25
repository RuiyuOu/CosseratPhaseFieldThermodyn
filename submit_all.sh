#!/bin/bash
# Submits a CHAIN of N restart-safe jobs per variant (A and F run in parallel), then post-processing after both chains.
# Each link starts only after the previous one ends for ANY reason (finished, 12 h time limit, crash) and resumes from the last
# checkpoint (saved every 0.5 ns = 50 steps). Links that find the run already complete exit immediately.
# Usage: bash submit_all.sh [N_links]   (default 4 => up to ~48 h of wall time per variant)
N=${1:-4}; mkdir -p logs
chain() { prev=$(sbatch --parsable run_dyn_$1.sbatch); for i in $(seq 2 $N); do prev=$(sbatch --parsable --dependency=afterany:$prev run_dyn_$1.sbatch); done; echo $prev; }
la=$(chain A); lf=$(chain F)
sbatch --dependency=afterany:$la:$lf run_post.sbatch
echo "submitted: A chain ends at job $la, F chain ends at job $lf; post depends on both"
