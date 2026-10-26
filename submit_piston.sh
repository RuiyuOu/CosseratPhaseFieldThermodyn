#!/bin/bash
# 6-level confined-piston ladder, 30 ns each, in parallel (array), plus N-1 follow-up arrays (afterany) for the 24 h limit.
# usage: bash submit_piston.sh [N_arrays=6]   optional env: TEND (ns, default 30), RAMP, HOLD (ns), VTKEVERY (steps)
N=${1:-6}; mkdir -p logs; export TEND RAMP HOLD VTKEVERY
prev=$(sbatch --parsable --export=ALL run_piston.sbatch)
for i in $(seq 2 $N); do prev=$(sbatch --parsable --export=ALL --dependency=afterany:$prev run_piston.sbatch); done
echo "piston ladder submitted; last array job $prev"
echo "outputs: out_pist_up<m/s>ms/pist_up<m/s>ms_s*.vtu (every 0.1 ns); a STALLED file marks a failed level"
