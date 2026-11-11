#!/bin/bash
# Clear the stop markers of crashed runs at one angle (a TSTOP above TBLOW = 5000 K is a blow-up and is rolled back), then resubmit.
# usage: ANGLE=0 bash resume.sh     env: MODE (mpi | serial; serial: DT 0.005), RUNTAG, DT, N (20), TBLOW, FORCE=1 (even if a run looks active)
set -uo pipefail
cd "$(dirname "$0")"
N=${N:-20}; TBLOW=${TBLOW:-5000}
export MODE=${MODE:-mpi}; [ "$MODE" = serial ] && export DT=${DT:-0.005}
export ANGLE=${ANGLE:-45}
UPMS=""; source ./tag.sh || exit 1
SPEEDS="25 50 100 150 250 400"
# never run two jobs on the same output directory: a run whose stage file was written in the last 20 min is still running
BUSY=""; for UP in $SPEEDS; do B=a${ANGLE}_up${UP}ms${SUFFIX}; [ -n "$(find out_$B/${B}_stage.txt -mmin -20 2>/dev/null)" ] && BUSY="$BUSY ${UP}"; done
if [ -n "$BUSY" ] && [ "${FORCE:-0}" != 1 ]; then
  echo "still running (stage file updated < 20 min ago) at$BUSY m/s for a${ANGLE}_up*ms${SUFFIX}: wait for those jobs or scancel them; FORCE=1 to override"; exit 1; fi
for UP in $SPEEDS; do
  B=a${ANGLE}_up${UP}ms${SUFFIX}; D=out_$B
  [ -d $D ] || { echo "$B: not started yet (will start fresh)"; continue; }
  if [ -e $D/${B}_TSTOP.txt ]; then
    TM=$(grep -o 'max T = [0-9.e+]*' $D/${B}_TSTOP.txt | awk '{print $4}')
    if awk -v t="${TM:-0}" -v c=$TBLOW 'BEGIN{exit !(t+0 > c)}'; then
      echo "$B: blow-up TSTOP (max T $TM K): rolling back"; bash recover_ck.sh $D || { echo "$B: recovery FAILED, not resubmitting this run"; continue; }
    else echo "$B: genuine TSTOP (max T $TM K <= $TBLOW K): left stopped"; continue; fi
  fi
  # NaN / non-converged stops save the clean state of the previous step: only the markers need to go
  for m in $D/${B}_NAN.txt $D/${B}_NOCONV.txt $D/STALLED; do [ -e "$m" ] && { mv "$m" "$m.old"; echo "$B: cleared $(basename $m)"; }; done
  CKG=$D/${B}_ck.txt; [ -s $CKG.manifest ] && CKG=$(sed -n 2p $CKG.manifest); [ -s $CKG ] && echo "$B: resumes from step $(head -1 $CKG)"
done
bash submit_speeds.sh $N
