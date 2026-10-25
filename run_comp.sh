#!/bin/bash
# 10 GPa uniaxial compression of variant A, started from the relaxed state of the finished run (skips the 100 relaxation steps).
# usage: bash run_comp.sh [amp_GPa=10] [tend_ns=10] [tscale=1] [tag=auto] [shear_GPa=1]   ; restart-safe like run_dyn.sh
# traction on the right edge = (-amp, +shear) GPa x envelope; shear=0 gives pure compression.
# tag=auto -> comp<AMP> (shear=0) or cs<AMP>_<SHEAR> (combined), so the cases never overwrite each other:
#   out_comp10/comp10_*.vtu   out_cs10_1/cs10_1_*.vtu   (plus the original shear case in out_A/A_dyn_*)
set -uo pipefail
AMP=${1:-10}; TEND=${2:-10}; TS=${3:-1}; SH=${5:-1}; TAG=${4:-auto}
if [ "$TAG" = auto ]; then if [ "$SH" = 0 ]; then TAG=comp${AMP}; else TAG=cs${AMP}_${SH}; fi; fi
D=out_$TAG; B=$TAG; mkdir -p $D
export FF_LOADPATH=${FF_LOADPATH:-$(dirname $(find $(dirname $(dirname $(which FreeFem++)))/lib* -name msh3.so 2>/dev/null | head -1) 2>/dev/null)}
done_p() { [ -s $D/${B}_log.csv ] && tail -1 $D/${B}_log.csv | cut -d, -f2 | awk -v e=$TEND '{exit !($1+0==$1 && $1>=e-0.011)}'; }
done_p && { echo "already complete"; exit 0; }
[ -e $D/STALLED ] && { echo "STALLED earlier (see $D/STALLED); fix the cause, then delete $D/STALLED to retry"; exit 1; }
maxstep() { [ -s $D/${B}_log.csv ] && awk -F, '$1+0==$1 && $1>m {m=$1} END{print m+0}' $D/${B}_log.csv || echo 0; }
S0=$(maxstep)
if [ ! -s $D/${B}_ck.txt ]; then
  [ -s out_A/A_dyn_relaxed_ck.txt ] || { echo "FAILED: out_A/A_dyn_relaxed_ck.txt missing (run variant A first)"; exit 1; }
  cp out_A/A_dyn_relaxed_ck.txt $D/${B}_ck.txt
fi
touch $D/$B.pvd
HDR="step,t_ns,load_GPa,iters,locres,KE_trans,KE_rot,E_elastic,E_stored,Q_heat,W_ext,D_cum,audit_residual,max_abs_slip,mean_abs_slip,max_dT,max_abs_dtheta_deg,mean_rhoS_um2,max_rhoG_um2,min_phi,neg_D_count"
[ -s $D/${B}_log.csv ] || echo "$HDR" > $D/${B}_log.csv   # solver writes the header only on a fresh (non-restart) start
DIRM=combo; [ "$SH" = 0 ] && DIRM=comp
ARGS="-variant A -j 0.01 -amp $AMP -dir $DIRM -shear $SH -tend $TEND -tscale $TS -dt 0.01 -muc 75.4 -nref 20 -vtkevery 5 -ckevery 50 -f0 0.005944 -kwca 0.10 -kwcs 1.05 -kwce 0.004 -kwcd 0.005"
FreeFem++ -nw -v 0 cosserat2d.edp $ARGS -out $D/$B -restart $D/${B}_ck.txt >> $D/$B.out 2>&1
RC=$?
S1=$(maxstep)
if done_p; then echo COMPLETE
elif [ "$S1" -le "$S0" ]; then
  { echo "STALLED: no new steps beyond step $S0 (solver exit code $RC) at $(date)"; echo "--- last lines of $D/$B.out ---"; tail -n 25 $D/$B.out; } | tee $D/STALLED
  exit 1
else echo "INCOMPLETE at step $S1 (exit code $RC): resubmit (restart-safe)"; fi
