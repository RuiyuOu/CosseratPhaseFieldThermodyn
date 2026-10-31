#!/bin/bash
# Roll a stopped run back to its newest rolling checkpoint so it can be resubmitted (e.g. after a spurious GB-spin blow-up:
# the TSTOP/NaN stop writes <tag>_ck.txt at the blown-up step, so resuming from it would start from the bad state).
# usage: bash recover_ck.sh out_<tag> [out_<tag2> ...]
#   - picks the newest <tag>_ckroll{0,1}.txt that is older than the step in <tag>_ck.txt
#   - copies it (and its .rhoM side file) to <tag>_ck.txt, keeping the bad one as <tag>_ck_bad<step>.txt
#   - moves the TSTOP / NAN / NOCONV / STALLED markers to *.old so the run scripts resume
#   log/diag/thermo rows past the checkpoint are trimmed by the run script itself on the next start.
set -uo pipefail
for D in "$@"; do
  D=${D%/}; B=${D#out_}
  [ -d "$D" ] || { echo "$D: no such directory"; continue; }
  CK=$D/${B}_ck.txt
  NBAD=$( [ -s "$CK" ] && head -1 "$CK" || echo 999999999 )
  best=""; nb=-1
  for r in $D/${B}_ckroll0.txt $D/${B}_ckroll1.txt; do
    [ -s "$r" ] || continue; n=$(head -1 "$r")
    if [ "$n" -lt "$NBAD" ] && [ "$n" -gt "$nb" ]; then best=$r; nb=$n; fi
  done
  [ -n "$best" ] || { echo "$D: no rolling checkpoint older than step $NBAD (was CKROLL on?)"; continue; }
  [ -s "$CK" ] && { mv "$CK" $D/${B}_ck_bad${NBAD}.txt; [ -f "$CK.rhoM" ] && mv "$CK.rhoM" $D/${B}_ck_bad${NBAD}.txt.rhoM; }
  cp "$best" "$CK"; [ -f "$best.rhoM" ] && cp "$best.rhoM" "$CK.rhoM"
  # drop the blown-up stop frame from the ParaView series (the frames after the checkpoint are rewritten on resume)
  [ -f $D/$B.pvd ] && { grep -v '_tstop[0-9]*\.vtu\|_stall[0-9]*\.vtu' $D/$B.pvd > $D/$B.pvd.tmp && mv $D/$B.pvd.tmp $D/$B.pvd; }
  for v in $D/${B}_tstop*.vtu $D/${B}_stall*.vtu; do [ -e "$v" ] && mv "$v" "$v.old"; done
  for m in $D/${B}_TSTOP.txt $D/${B}_NAN.txt $D/${B}_NOCONV.txt $D/STALLED; do [ -e "$m" ] && mv "$m" "$m.old"; done
  echo "$D: checkpoint rolled back from step $NBAD to step $nb ($(basename $best)); markers moved to *.old"
done
