#!/bin/bash
# Self-test (serial): AD element map vs central differences (-adcheck 1, every element), then the global tangent vs central differences
# (-tancheck 2) and the energy derivatives (-check 1) over a few loaded steps from the relaxation checkpoint of the set-up.
# usage: bash checks/run_tangent_check.sh [steps=3]     env: SPEED (400), ANGLE (0), TCSTEP, OLD (0)
#   The v2 solver with direct assembly, the exact slip law and the acceptance tests of production runs, and the exact KWC Hessian
#   (-kwctan 0: the derivative of the primal residual; -kwctan 2 is not, away from convergence). OLD=0 (default): production parameters, multiplier (held fixed by the FD check), relax_checkpoint_v2.txt;
#   OLD=1: rotational inertia (j = 0.01), otherwise the same, relax_checkpoint_j.txt
# measured (coarse mesh): adcheck |FD - AD| ~1e-7 of max |AD| (the FD error), values ~1e-13; tangent ~1e-7 (h^2 down to ~1e-8), KWC tangent
#           ~1e-11, multiplier blocks C, B, E at the FD error; tau = M:B to the FD error; internal force = dE/dU to ~6 digits. The
#           finite-increment driving-force comparison is informational (Gauss-point energy, finite exponential update).
# The last line is decided by checks/eval_checks.py: SELF-TEST PASSED / FAILED (the script exits with 1 on a failure): both FreeFEM
# exit statuses 0, the loaded run complete (final step NR + NS, no _NOCONV/_NAN/_TSTOP marker, no ERROR line), every record present
# exactly once, finite and within its limit. checks/test_eval_checks.py (run first) tests the evaluator on fixtures.
set -euo pipefail
cd "$(dirname "$0")/.."
NS=${1:-3}; SPEED=${SPEED:-400}; export ANGLE=${ANGLE:-0}
export OLD=${OLD:-0}; SRC=cosserat_kwc_v2.edp; CHK=cosserat_kwc_v2_check.edp
if [ "$OLD" = 1 ]; then RCK=relax_checkpoint_j.txt; else RCK=relax_checkpoint_v2.txt; fi; XO="-kwctan 0 -almc 4 -asm 1"
XO="$XO -slipexact 1 -localtol 1e-10 -betamax 0.785398"   # the production slip law and acceptance tests (as run.sh passes them)
[ -s $RCK ] || { echo "$RCK missing: run '$([ "$OLD" = 1 ] && echo "OLD=1 ")MODE=serial bash run_relax.sh' first"; exit 1; }
python3 checks/test_eval_checks.py > /dev/null || { echo "FAILED: the pass/fail evaluator gives a wrong verdict on its fixtures (python3 checks/test_eval_checks.py)"; exit 1; }
bash plugin/build_plugin.sh
# generated in the repository root: FreeFem++ reads '-cd' (a solver option in PHYS) as 'change to the script's directory', so the
# script must sit next to mesh.msh, the grain files and kwc2d.so
python3 checks/make_checks.py $SRC $CHK
MODE=serial; UPMS=$SPEED; source ./tag.sh
D=checks/out_tancheck$([ "$OLD" = 1 ] && echo _j || echo _v2); rm -rf $D; mkdir -p $D; cp $RCK $D/t_ck.txt
NR=$(head -1 $RCK); NRX=${NRELAX:-100}; [ "$NR" -gt "$NRX" ] && NRX=$NR   # loaded steps start after max(checkpoint step, -nrelax)
# the exit status of each FreeFEM run is kept (set +e around them) and is part of the pass/fail decision
set +e
FreeFem++ -nw -v 0 $CHK $PHYS $XO -adcheck 1 -nrelax 0 -out $D/ad > $D/ad.out 2>&1; RCAD=$?
FreeFem++ -nw -v 0 $CHK -pistonup $UN -pistonut $UT -tend $(awk -v n=$NS 'BEGIN{print n*0.01}') -vtkevery 0 $PHYS $XO \
  -tancheck 2 -tancheckstep ${TCSTEP:-$((NRX + NS))} -check 1 -out $D/t -restart $D/t_ck.txt > $D/t.out 2>&1; RCRUN=$?
set -e
grep -E "adcheck|ERROR|error" $D/ad.out || true
grep -E "TANGENT CHECK|FD check|done\.|ERROR|error" $D/t.out || true
echo "exit status: adcheck $RCAD, loaded run $RCRUN"
# pass/fail (checks/eval_checks.py): statuses, completion (final step NR + NS, no failure marker), every record present once, finite, within its limit
python3 checks/eval_checks.py --ad $D/ad.out --ad-status $RCAD --run $D/t.out --run-status $RCRUN --log $D/t_log.csv \
  --expect-step $((NRX + NS)) --prefix $D/t --multiplier 1
