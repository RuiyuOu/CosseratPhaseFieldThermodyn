#!/bin/bash
# Roll a stopped run back to its newest rolling checkpoint so it can be resubmitted (e.g. after a spurious GB-spin blow-up:
# the TSTOP/NaN stop writes <tag>_ck.txt at the blown-up step, so resuming from it would start from the bad state).
# usage: bash recover_ck.sh out_<tag> [out_<tag2> ...]
#   - picks the newest <tag>_ckroll{0,1}.txt that is older than the step in <tag>_ck.txt
#   - makes it the committed checkpoint <tag>_ck.txt (checkpoints are generations selected by <tag>_ck.txt.manifest: the rolling checkpoint's
#     generation is hard-linked as a new generation of <tag>_ck.txt and the manifest is replaced; the bad one's manifest is archived by copy as <tag>_ck_bad<step>.txt.manifest)
#   - moves the TSTOP / NAN / NOCONV / STALLED markers to *.old so the run scripts resume
#   log/diag/thermo rows past the checkpoint are trimmed by the run script itself on the next start.
set -uo pipefail
rc=0
for D in "$@"; do
  D=${D%/}; B=${D#out_}
  [ -d "$D" ] || { echo "$D: no such directory"; continue; }
  CK=$D/${B}_ck.txt
  # committed generation of a checkpoint name: the manifest's generation if there is a manifest, else the plain file
  gen_of() { if [ -s "$1.manifest" ]; then sed -n 2p "$1.manifest"; else echo "$1"; fi; }
  GCK=$(gen_of "$CK")
  NBAD=$( [ -s "$GCK" ] && head -1 "$GCK" || echo 999999999 )
  best=""; nb=-1
  for r in $D/${B}_ckroll0.txt $D/${B}_ckroll1.txt; do
    g=$(gen_of "$r"); [ -s "$g" ] || continue; n=$(head -1 "$g")
    if [ "$n" -lt "$NBAD" ] && [ "$n" -gt "$nb" ]; then best=$r; nb=$n; fi
  done
  [ -n "$best" ] || { echo "$D: no rolling checkpoint older than step $NBAD (was CKROLL on?)" >&2; rc=1; continue; }
  G=$(gen_of "$best"); NG=$CK.g${nb}_rb_$(date +%s)_$$
  # 1. prepare a uniquely named generation; the live manifest is not touched until every step has succeeded
  fail=0
  for suf in "" .rhoM .gnd .lam .thk .plate; do
    if [ -e "$G$suf" ]; then ln -f "$G$suf" "$NG$suf" 2>/dev/null || cp -f "$G$suf" "$NG$suf" || { fail=1; break; }; fi
  done
  # required components: the state and its rho_M file must be present and non-empty, and the step must match
  if [ $fail = 0 ]; then
    { [ -s "$NG" ] && [ -s "$NG.rhoM" ] && [ "$(head -1 "$NG")" = "$nb" ]; } || fail=1
    for suf in .gnd .lam .thk .plate; do [ -e "$G$suf" ] && { [ -s "$NG$suf" ] || fail=1; }; done
  fi
  GC=0; [ -s "$best.manifest" ] && GC=$(sed -n 4p "$best.manifest")
  if [ $fail = 0 ]; then
    printf '%s\n%s\n-\n%s\n' "$nb" "$NG" "${GC:-0}" > "$CK.manifest.tmp" && [ -s "$CK.manifest.tmp" ] || fail=1
  fi
  # 2. archive the old pointer by copy, then commit by one rename
  if [ $fail = 0 ] && [ -s "$CK.manifest" ]; then cp -f "$CK.manifest" "$D/${B}_ck_bad${NBAD}.txt.manifest" || fail=1; fi
  if [ $fail = 0 ]; then mv -f "$CK.manifest.tmp" "$CK.manifest" || fail=1; fi
  if [ $fail != 0 ]; then
    rm -f "$NG" "$NG".rhoM "$NG".gnd "$NG".lam "$NG".thk "$NG".plate "$CK.manifest.tmp"
    echo "$D: RECOVERY FAILED (file operation); live manifest and stop markers left unchanged" >&2
    rc=1; continue
  fi
  # 3. committed: refresh the convenience links (best effort; the manifest is authoritative) and clear the markers
  for suf in "" .rhoM .gnd .lam .thk .plate; do
    if [ -e "$NG$suf" ]; then ln -f "$NG$suf" "$CK$suf" 2>/dev/null || echo "$D: warning: could not refresh $CK$suf (manifest is authoritative)" >&2; else rm -f "$CK$suf"; fi
  done
  # drop the blown-up stop frame from the ParaView series (the frames after the checkpoint are rewritten on resume)
  [ -f $D/$B.pvd ] && { grep -v '_tstop[0-9]*\.vtu\|_stall[0-9]*\.vtu' $D/$B.pvd > $D/$B.pvd.tmp && mv $D/$B.pvd.tmp $D/$B.pvd; }
  for v in $D/${B}_tstop*.vtu $D/${B}_stall*.vtu; do [ -e "$v" ] && mv "$v" "$v.old"; done
  for m in $D/${B}_TSTOP.txt $D/${B}_NAN.txt $D/${B}_NOCONV.txt $D/STALLED; do [ -e "$m" ] && mv "$m" "$m.old"; done
  echo "$D: checkpoint rolled back from step $NBAD to step $nb ($(basename $best)); markers moved to *.old"
done
exit $rc
