#!/bin/bash
# Resume the whole serial variant-A piston ladder (25 50 100 150 250 400 m/s) after crashes: clears the stop markers of
# crashed levels (rolling a blown-up TSTOP back to its last rolling checkpoint), then submits the chain again
# (run_all_A.sh: N chained 24 h arrays + post-processing). Finished levels exit at once.
# usage (from hpc_dynamic/):
#   RUNTAG=gm bash resume_all_A.sh                       # DT = 0.01 ns, NREF = 20   -> out_pistM3v2RcV_gm_up<speed>ms
#   RUNTAG=gm DT=0.005 NREF=10 bash resume_all_A.sh      # half dt, NREF 10           -> out_pistM3v2RcV_dt0.005_nr10_gm_up<speed>ms
#   MPI=1: the MPI runs (tag suffix _mpi, submits run_all_A_mpi.sh; NP, MEM as there); GNDPHI2=0: the runs without the phi^2 GND weight (tag G)
#   optional: N=20 chained arrays (default 20); TBLOW=5000: a TSTOP above this max T (K) counts as a numerical blow-up and is rolled
#   back (below it the run genuinely reached -tstop, i.e. melting, and is left stopped); FORCE=1 submits even if a level of this tag is still running.
set -uo pipefail
cd "$(dirname "$0")"
export DT=${DT:-0.01} NREF=${NREF:-20}
N=${N:-20}; TBLOW=${TBLOW:-5000}
G=""; [ "${GNDPHI2:-1}" = 0 ] && G=G; export GNDPHI2=${GNDPHI2:-1}
SFX=M3v2${G}RcV; { [ "$DT" != 0.01 ] || [ "$NREF" != 20 ]; } && SFX=${SFX}_dt${DT}_nr${NREF}; [ -n "${RUNTAG:-}" ] && SFX=${SFX}_${RUNTAG}
[ "${MPI:-0}" = 1 ] && [ "${SHARETAG:-0}" != 1 ] && SFX=${SFX}_mpi
# never run two jobs on the same output directory: a level whose stage file was written in the last 20 min is still running
BUSY=""; for UP in 25 50 100 150 250 400; do f=out_pist${SFX}_up${UP}ms/pist${SFX}_up${UP}ms_stage.txt; [ -n "$(find $f -mmin -20 2>/dev/null)" ] && BUSY="$BUSY ${UP}"; done
if [ -n "$BUSY" ] && [ "${FORCE:-0}" != 1 ]; then
  echo "still running (stage file updated < 20 min ago) at u_p =$BUSY m/s for tag pist${SFX}: wait for those jobs or scancel them; FORCE=1 to override"; exit 1; fi
for UP in 25 50 100 150 250 400; do
  B=pist${SFX}_up${UP}ms; D=out_$B
  [ -d $D ] || { echo "$B: not started yet (will start fresh)"; continue; }
  if [ -e $D/${B}_TSTOP.txt ]; then
    TM=$(grep -o 'max T = [0-9.e+]*' $D/${B}_TSTOP.txt | awk '{print $4}')
    if awk -v t="${TM:-0}" -v c=$TBLOW 'BEGIN{exit !(t+0 > c)}'; then
      echo "$B: blow-up TSTOP (max T $TM K): rolling back"; bash recover_ck.sh $D
    else echo "$B: genuine TSTOP (max T $TM K <= $TBLOW K): left stopped"; continue; fi
  fi
  # NaN / non-converged stops save the clean state of the previous step: only the markers need to go
  for m in $D/${B}_NAN.txt $D/${B}_NOCONV.txt $D/STALLED; do [ -e "$m" ] && { mv "$m" "$m.old"; echo "$B: cleared $(basename $m)"; }; done
  [ -s $D/${B}_ck.txt ] && echo "$B: resumes from step $(head -1 $D/${B}_ck.txt)"
done
RS=run_piston_mob.sh; [ "${MPI:-0}" = 1 ] && RS=run_piston_mob_mpi.sh; grep -q spinbranch $RS || echo "WARNING: $RS has no -spinbranch (apply the spin fix first, or the same blow-up repeats)"
if [ "${MPI:-0}" = 1 ]; then bash run_all_A_mpi.sh $N; else bash run_all_A.sh $N; fi
