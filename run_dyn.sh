#!/bin/bash
# usage: bash run_dyn.sh <A|F> ; restart-safe: resumes from <out>_ck.txt, exits immediately if the run already reached 15 ns
set -uo pipefail
V=${1:?variant A or F}; D=out_$V; B=${V}_dyn; mkdir -p $D
export FF_LOADPATH=${FF_LOADPATH:-$(dirname $(find $(dirname $(dirname $(which FreeFem++)))/lib* -name msh3.so 2>/dev/null | head -1) 2>/dev/null)}
done_p() { [ -s $D/${B}_log.csv ] && tail -1 $D/${B}_log.csv | cut -d, -f2 | awk '{exit !($1+0==$1 && $1>=14.99)}'; }
done_p && { echo "already complete"; exit 0; }
[ -e $D/${B}_TSTOP.txt ] && { echo "ended at the temperature limit:"; cat $D/${B}_TSTOP.txt; exit 0; }
[ -e $D/STALLED ] && { echo "STALLED earlier (see $D/STALLED); fix the cause, then delete $D/STALLED to retry"; exit 1; }
maxstep() { [ -s $D/${B}_log.csv ] && awk -F, '$1+0==$1 && $1>m {m=$1} END{print m+0}' $D/${B}_log.csv || echo 0; }
S0=$(maxstep)
ARGS="-j 0.01 -amp 0.06 -dir shear -tend 15 -dt 0.01 -muc 75.4 -nref 20 -vtkevery 5 -ckevery 50 -f0 0.005944 -kwca 0.10 -kwcs 1.05 -kwce 0.004 -kwcd 0.005"
touch $D/$B.pvd; R=""; [ -s $D/${B}_ck.txt ] && R="-restart $D/${B}_ck.txt"
FreeFem++ -nw -v 0 cosserat2d.edp -variant $V $ARGS -out $D/$B $R >> $D/$B.out 2>&1
done_p && echo "COMPLETE" || echo "INCOMPLETE: resubmit (restart-safe)"
