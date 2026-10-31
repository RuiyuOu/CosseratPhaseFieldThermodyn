#!/bin/bash
# Variant F only: 6-level piston ladder (30 ns), in parallel, plus N-1 follow-up arrays (afterany) for the 24 h limit.
# usage: bash submit_piston_F.sh [N_arrays=3]    optional env: TEND, RAMP, HOLD, VTKEVERY, TSTOP
N=${1:-3}; mkdir -p logs; export TEND RAMP HOLD VTKEVERY TSTOP
prev=$(sbatch --parsable --export=ALL run_piston_F.sbatch)
for i in $(seq 2 $N); do prev=$(sbatch --parsable --export=ALL --dependency=afterany:$prev run_piston_F.sbatch); done
echo "variant F piston ladder submitted; last array job $prev"
echo "outputs: out_pistF_up<m/s>ms/ ; post-process with: sbatch run_post_piston_F.sbatch"
