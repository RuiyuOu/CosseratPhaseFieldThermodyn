#!/bin/bash
# Self-test (serial): tangent vs central differences (-tancheck 2) and energy derivatives (-check 1), a few loaded steps from relax_checkpoint.txt.
# usage: bash checks/run_tangent_check.sh [steps=3]     env: SPEED (400), ANGLE (0), TCSTEP
# expected: tangent error ~1e-7 (h^2 down to ~1e-8), KWC tangent ~1e-11, FD checks agreeing to ~1e-6
set -euo pipefail
cd "$(dirname "$0")/.."
NS=${1:-3}; SPEED=${SPEED:-400}; export ANGLE=${ANGLE:-0}
[ -s relax_checkpoint.txt ] || { echo "relax_checkpoint.txt missing: run 'MODE=serial bash run_relax.sh' first"; exit 1; }
python3 checks/make_checks.py cosserat_kwc.edp checks/cosserat_kwc_check.edp
MODE=serial; UPMS=$SPEED; source ./tag.sh
D=checks/out_tancheck; rm -rf $D; mkdir -p $D; cp relax_checkpoint.txt $D/t_ck.txt
NR=$(head -1 relax_checkpoint.txt)
FreeFem++ -nw -v 0 checks/cosserat_kwc_check.edp -pistonup $UN -pistonut $UT -tend $(awk -v n=$NS 'BEGIN{print n*0.01}') -vtkevery 0 $PHYS \
  -tancheck 2 -tancheckstep ${TCSTEP:-$((NR + NS))} -check 1 -out $D/t -restart $D/t_ck.txt 2>&1 | grep -E "TANGENT CHECK|FD check|done\.|ERROR|error"
