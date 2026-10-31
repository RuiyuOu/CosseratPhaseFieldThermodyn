#!/bin/bash
# GND/recovery piston ladder (variant A x 6 speeds; defaults GNDPHI2=0, RECOV=2 -> tag pistGR), 30 ns each, plus N-1 follow-up arrays (afterany) for the 24 h limit.
# usage: bash submit_piston_gnd.sh [N_arrays=6]    optional env: TEND, RAMP, HOLD, VTKEVERY, TSTOP, GNDPHI2 (default 0), RECOV (default 2), CA (default sqrt 10 um), CD (default 100), DIAG (default 1)
# needs out_A/A_dyn_relaxed_ck.txt (same relaxed start as the original A ladder)
N=${1:-6}; mkdir -p logs; export TEND RAMP HOLD VTKEVERY TSTOP GNDPHI2 RECOV CA CD DIAG
[ -s out_A/A_dyn_relaxed_ck.txt ] || { echo "missing out_A/A_dyn_relaxed_ck.txt"; exit 1; }
prev=$(sbatch --parsable --export=ALL run_piston_gnd.sbatch)
for i in $(seq 2 $N); do prev=$(sbatch --parsable --export=ALL --dependency=afterany:$prev run_piston_gnd.sbatch); done
echo "GND/recovery piston ladder submitted; last array job $prev"
echo "outputs: out_pist<G|R|GR>_up<m/s>ms/ (default GR)  (frames every 0.1 ns, <tag>_diag.csv = worst element per step); a STALLED file marks a failed level"
