#!/bin/bash
# Confined piston (cold-spray analogue), variant A, resumes from out_A/A_dyn_relaxed_ck.txt. Restart-safe, with the STALLED guard.
# usage: bash run_piston.sh <up_m_per_s> [tend_ns=30] [ramp_ns=0.5] [hold_ns=2.0] [vtkevery=10]
#   piston velocity u_p (m/s) is applied to the right edge (u_x = -d(t)); wall u_x=0 on the left; u_y=0 on top and bottom.
#   total piston travel = u_p*(ramp+hold); frames every vtkevery*dt = 0.1 ns by default (dt = 0.01 ns).
set -uo pipefail
UPMS=${1:?piston velocity in m/s}; TEND=${2:-30}; TR=${3:-0.5}; TH=${4:-2.0}; VE=${5:-10}
UP=$(awk -v u=$UPMS 'BEGIN{printf "%g", u/1000}')           # um/ns
TAG=pist_up${UPMS}ms; D=out_$TAG; B=$TAG; mkdir -p $D
export FF_LOADPATH=${FF_LOADPATH:-$(dirname $(find $(dirname $(dirname $(which FreeFem++)))/lib* -name msh3.so 2>/dev/null | head -1) 2>/dev/null)}
done_p() { [ -s $D/${B}_log.csv ] && tail -1 $D/${B}_log.csv | cut -d, -f2 | awk -v e=$TEND '{exit !($1+0==$1 && $1>=e-0.011)}'; }
done_p && { echo "already complete"; exit 0; }
[ -e $D/${B}_TSTOP.txt ] && { echo "ended at the temperature limit:"; cat $D/${B}_TSTOP.txt; exit 0; }
[ -e $D/STALLED ] && { echo "STALLED earlier (see $D/STALLED); fix the cause, then delete $D/STALLED to retry"; exit 1; }
maxstep() { [ -s $D/${B}_log.csv ] && awk -F, '$1+0==$1 && $1>m {m=$1} END{print m+0}' $D/${B}_log.csv || echo 0; }
S0=$(maxstep)
if [ ! -s $D/${B}_ck.txt ]; then
  [ -s out_A/A_dyn_relaxed_ck.txt ] || { echo "FAILED: out_A/A_dyn_relaxed_ck.txt missing (run variant A first)"; exit 1; }
  cp out_A/A_dyn_relaxed_ck.txt $D/${B}_ck.txt
fi
touch $D/$B.pvd
HDR="step,t_ns,load_GPa,iters,locres,KE_trans,KE_rot,E_elastic,E_stored,Q_heat,W_ext,D_cum,audit_residual,max_abs_slip,mean_abs_slip,max_dT,max_abs_dtheta_deg,mean_rhoS_um2,max_rhoG_um2,min_phi,neg_D_count"
[ -s $D/${B}_log.csv ] || echo "$HDR" > $D/${B}_log.csv
ARGS="-variant A -j 0.01 -amp 1 -dir piston -pistonup $UP -ptr $TR -pth $TH -tend $TEND -dt 0.01 -muc 75.4 -nref 20 -vtkevery $VE -ckevery 50 -tstop ${TSTOP:-1400} -f0 0.005944 -kwca 0.10 -kwcs 1.05 -kwce 0.004 -kwcd 0.005"
FreeFem++ -nw -v 0 cosserat2d.edp $ARGS -out $D/$B -restart $D/${B}_ck.txt >> $D/$B.out 2>&1
RC=$?
S1=$(maxstep)
if done_p; then echo COMPLETE
elif [ -e $D/${B}_TSTOP.txt ]; then echo "ENDED at the temperature limit:"; cat $D/${B}_TSTOP.txt
elif [ "$S1" -le "$S0" ]; then
  { echo "STALLED: no new steps beyond step $S0 (solver exit code $RC) at $(date)"; echo "--- last lines of $D/$B.out ---"; tail -n 25 $D/$B.out; } | tee $D/STALLED
  exit 1
else echo "INCOMPLETE at step $S1 (exit code $RC): resubmit (restart-safe)"; fi
