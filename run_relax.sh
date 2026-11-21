#!/bin/bash
# Relaxation stage: makes the relaxation checkpoint of the set-up (relax_checkpoint_v2.txt, OLD=1: relax_checkpoint_j.txt) if it is missing, does nothing if it exists with the same parameters, refuses otherwise.
# usage: bash run_relax.sh     env: MODE (mpi | serial), NP, RELAXCK, RELAXCHECK, physics variables of tag.sh (must match the runs)
set -uo pipefail
cd "$(dirname "$0")"
RELAXONLY=1 bash ./run.sh 0 0 ${RAMP:-0.5} ${HOLD:-2.0} ${VTKEVERY:-10}
