#!/bin/bash
# chain of N restart-safe jobs for the compression test: bash submit_comp.sh [N=6]
# optional env: AMP (GPa), SHEAR (GPa, default 1), TEND (ns), TSCALE (load time scale), TAG (output tag) e.g.  AMP=10 TSCALE=0.5 TAG=C10fast bash submit_comp.sh
N=${1:-6}; mkdir -p logs; export AMP TEND TSCALE TAG SHEAR
prev=$(sbatch --parsable --export=ALL run_comp.sbatch); for i in $(seq 2 $N); do prev=$(sbatch --parsable --export=ALL --dependency=afterany:$prev run_comp.sbatch); done; echo "chain ends at job $prev"
