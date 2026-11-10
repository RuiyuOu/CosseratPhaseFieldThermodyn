#!/bin/bash
# Resume the whole serial variant-A piston ladder (25 50 100 150 250 400 m/s) after crashes: clears the stop markers of
# crashed levels (rolling a blown-up TSTOP back to its last rolling checkpoint), then submits the chain again
# (run_all_A.sh: N chained 24 h arrays + post-processing). Finished levels exit at once.
# usage (from hpc_dynamic/):
#   RUNTAG=gm bash resume_all_A.sh                       # DT = 0.01 ns, NREF = 20   -> out_pistM3v2GRcV_gc_ne_pm1e-2_es1_gm_a45_up<speed>ms
#   RUNTAG=gm DT=0.005 NREF=10 bash resume_all_A.sh      # half dt, NREF 10           -> out_pistM3v2GRcV_dt0.005_nr10_gc_ne_pm1e-2_es1_gm_a45_up<speed>ms
#   MPI=1: the MPI runs (tag suffix _mpi, submits run_all_A_mpi.sh; NP, MEM as there); GND / legacy GNDCURL, GNDPHI2, KWCEDIFF, ETAPHIMIN, ETASTAR, CODE, EOS ...: as given to the run scripts (all part of the tag)
#   SYMFOLD / GNDCURL / SPINMAX as given to the run scripts (they are part of the tag)
#   ANGLE (default 45, as in the run scripts; ANGLE=confined for the old confined-piston runs, e.g. the bk series)
#   optional: N=20 chained arrays (default 20); TBLOW=5000: a TSTOP above this max T (K) counts as a numerical blow-up and is rolled
#   back (below it the run genuinely reached -tstop, i.e. melting, and is left stopped); FORCE=1 submits even if a level of this tag is still running.
set -uo pipefail
cd "$(dirname "$0")"
export DT=${DT:-0.01} NREF=${NREF:-20}
N=${N:-20}; TBLOW=${TBLOW:-5000}
if [ -n "${GNDCURL:-}" ] || [ -n "${GNDPHI2:-}" ]; then GP=${GNDPHI2:-1}; GCL=${GNDCURL:-0}   # legacy flags win
else GND=${GND:-curl}; case $GND in curl) GP=0; GCL=1;; grad) GP=0; GCL=0;; phigrad) GP=1; GCL=0;; *) echo "GND must be curl|grad|phigrad" >&2; exit 1;; esac; fi
SFX="M${DRAGRHO:-3}v${VC:-2}"; [ "$GP" = 0 ] && SFX=${SFX}G; [ "${RECOV:-2}" != 0 ] && SFX=${SFX}R; { [ "${DISSIP:-1}" != 0 ] || [ "${CDSTORE:-1}" != 0 ]; } && SFX=${SFX}c; case ${EOS:-vinet} in vinet) SFX=${SFX}V;; murnaghan) SFX=${SFX}Mu;; esac; { [ "${ELINT:-$([ "${CODE:-stable}" = stable ] && echo gauss || echo mean)}" = gauss ] && SFX=${SFX}Q${NQ:-3}; }; { [ "${DT:-0.01}" != 0.01 ] || [ "${NREF:-20}" != 20 ]; } && SFX=${SFX}_dt${DT:-0.01}_nr${NREF:-20}; [ "${SPINMAX:-0}" != 0 ] && SFX=${SFX}_sm${SPINMAX}; [ "${SYMFOLD:-0}" != 0 ] && SFX=${SFX}_sf${SYMFOLD}; [ "$GCL" != 0 ] && SFX=${SFX}_gc; [ "${CODE:-stable}" = exp ] && SFX=${SFX}_exp; [ "${MUCDUAL:-0}" != 0 ] && SFX=${SFX}_du; [ "${MUCDUAL:-0}" != 0 ] && [ "${DUALR:-0}" != 0 ] && SFX=${SFX}_dr${DUALR}; [ "${KWCEDIFF:-0}" = 0 ] && SFX=${SFX}_ne; [ "${KWCE:-0.004}" != 0.004 ] && SFX=${SFX}_e${KWCE}; [ "${KWCREG:-0}" != 0 ] && SFX=${SFX}_kg; [ "${KWCREG:-0}" != 0 ] && [ "${KWCGAM:-0}" != 0 ] && SFX=${SFX}${KWCGAM}; [ "${KWCFIELD:-0}" != 0 ] && SFX=${SFX}_kf; [ "${KWCONEWAY:-0}" != 0 ] && SFX=${SFX}_ow; PMN=${ETAPHIMIN:-1e-2}; [ "$PMN" != 0 ] && SFX=${SFX}_pm${PMN}; [ "${QPHI:-100e3}" != 100e3 ] && SFX=${SFX}_q${QPHI}; [ "${ETAPHI:-0.1}" != 0.1 ] && SFX=${SFX}_ep${ETAPHI}; EST=${ETASTAR:-1}; [ "$EST" != 1e-3 ] && SFX=${SFX}_es${EST}   # same tag logic as run_post_piston_mob.sbatch
[ -n "${RUNTAG:-}" ] && SFX=${SFX}_${RUNTAG}
[ "${MPI:-0}" = 1 ] && [ "${SHARETAG:-0}" != 1 ] && SFX=${SFX}_mpi
ATAG=""; [ "${ANGLE:-45}" != confined ] && ATAG=_a${ANGLE:-45}; export ANGLE=${ANGLE:-45}   # ANGLE=confined: the old confined-piston runs
# never run two jobs on the same output directory: a level whose stage file was written in the last 20 min is still running
BUSY=""; for UP in 25 50 100 150 250 400; do f=out_pist${SFX}${ATAG}_up${UP}ms/pist${SFX}${ATAG}_up${UP}ms_stage.txt; [ -n "$(find $f -mmin -20 2>/dev/null)" ] && BUSY="$BUSY ${UP}"; done
if [ -n "$BUSY" ] && [ "${FORCE:-0}" != 1 ]; then
  echo "still running (stage file updated < 20 min ago) at u_p =$BUSY m/s for tag pist${SFX}${ATAG}: wait for those jobs or scancel them; FORCE=1 to override"; exit 1; fi
for UP in 25 50 100 150 250 400; do
  B=pist${SFX}${ATAG}_up${UP}ms; D=out_$B
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
RS=run_piston_mob.sh; [ "${MPI:-0}" = 1 ] && RS=run_piston_mob_mpi.sh; grep -q spinbranch $RS || echo "WARNING: $RS has no -spinbranch (apply the spin fix first, or the same blow-up repeats)"
if [ "${MPI:-0}" = 1 ]; then bash run_all_A_mpi.sh $N; else bash run_all_A.sh $N; fi
