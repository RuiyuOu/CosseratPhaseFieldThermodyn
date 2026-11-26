#!/bin/bash
# push_run_data.sh -- commit one run's output to a data branch of this repository and push it, in parts.
#
#   bash push_run_data.sh <tag> [all|meta|vtu] [first VTU group]
#     all  (default): first every file except the VTUs, then the VTUs in groups of STEPS loaded steps up to KMAX
#     meta           : only the files except the VTUs
#     vtu  [g]       : only the VTU groups, starting at group g (0 = steps 0-499 and the relaxed frame)
#   env: BRANCH (data)  REMOTE (origin)  STEPS (500)  KMAX (3000 = 30 ns at dt 0.01)  MAXMB (95: larger files are skipped, GitHub
#        rejects files over 100 MB)
#
# What goes in: out_<tag>/ (every file; VTUs by group), deliver_<tag>/grains/ and deliver_<tag>/tau/ (not the GIFs and the data zip, which
# repeat the run folder). Each part is one commit and one push, so a failed push loses nothing: rerun with the group it stopped at.
#
# The commits go to the branch BRANCH (default 'data'), built with a separate index: the branch you have checked out, its index and your
# working tree are not touched, .gitignore (which excludes out_*/ and deliver_*/) does not apply, and the code history stays separate.
# Fetch it elsewhere with  git fetch origin data:data  (or a partial clone: git clone --filter=blob:none --no-checkout -b data <url>).
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
TAG=${1:?usage: bash push_run_data.sh <tag> [all|meta|vtu] [first VTU group]}
WHAT=${2:-all}; G0=${3:-0}
BR=${BRANCH:-data}; REMOTE=${REMOTE:-origin}; STEPS=${STEPS:-500}; KMAX=${KMAX:-3000}; MAXMB=${MAXMB:-95}
D=out_$TAG
case "$WHAT" in all|meta|vtu) ;; *) echo "second argument: all, meta or vtu"; exit 1;; esac
[ -d "$D" ] || { echo "no folder $D"; exit 1; }
git rev-parse --git-dir >/dev/null 2>&1 || { echo "not inside the git repository"; exit 1; }
CUR=$(git symbolic-ref --short -q HEAD || true)
[ "$BR" != "$CUR" ] || { echo "BRANCH=$BR is the branch checked out here: use another branch (default 'data') so the working tree is not disturbed"; exit 1; }
git remote get-url "$REMOTE" >/dev/null 2>&1 || { echo "no remote '$REMOTE'"; exit 1; }

export GIT_INDEX_FILE="$(git rev-parse --absolute-git-dir)/index.rundata.$$"
trap 'rm -f "$GIT_INDEX_FILE"' EXIT
# start from the tip of the data branch: local, else the remote one, else a new branch with no history
if git rev-parse -q --verify "refs/heads/$BR" >/dev/null; then P=$(git rev-parse "refs/heads/$BR")
elif git ls-remote --exit-code --heads "$REMOTE" "$BR" >/dev/null 2>&1; then git fetch -q "$REMOTE" "$BR:refs/heads/$BR"; P=$(git rev-parse "refs/heads/$BR")
else P=""; fi
if [ -n "$P" ]; then git read-tree "$P"; echo "branch $BR at ${P:0:12}"; else git read-tree --empty; echo "branch $BR: new (no history)"; fi

LIST=$(mktemp); trap 'rm -f "$GIT_INDEX_FILE" "$LIST" "$LIST.ok"' EXIT
part() {   # $1 = commit message; files in $LIST
  : > "$LIST.ok"
  while IFS= read -r f; do
    [ -f "$f" ] || continue
    if [ "$(stat -c %s "$f")" -gt $((MAXMB*1024*1024)) ]; then echo "  SKIPPED (> $MAXMB MB): $f"; else echo "$f" >> "$LIST.ok"; fi
  done < "$LIST"
  local n; n=$(wc -l < "$LIST.ok")
  [ "$n" -gt 0 ] || { echo "  no files"; return 0; }
  echo "  $n files, $(xargs -a "$LIST.ok" -d '\n' du -ch | tail -1 | cut -f1)"
  xargs -a "$LIST.ok" -d '\n' git update-index --add --
  local T C; T=$(git write-tree)
  if [ -n "$P" ] && [ "$T" = "$(git rev-parse "$P^{tree}")" ]; then echo "  already on $BR: nothing to commit"; else
    C=$(git commit-tree "$T" ${P:+-p "$P"} -m "$1")
    git update-ref "refs/heads/$BR" "$C" "$P"
    P=$C; echo "  commit ${C:0:12}"
  fi
  git -c http.postBuffer=524288000 push -q "$REMOTE" "refs/heads/$BR:refs/heads/$BR" && echo "  pushed to $REMOTE/$BR"
}

if [ "$WHAT" != vtu ]; then
  echo "[$TAG] files except the VTUs"
  find "$D" "deliver_$TAG/grains" "deliver_$TAG/tau" -type f ! -name '*.vtu' ! -name '*.tmp' ! -name '*.lock' 2>/dev/null | sort > "$LIST"
  part "$TAG: run files except the VTUs"
fi
if [ "$WHAT" != meta ]; then
  NG=$(( (KMAX + STEPS - 1)/STEPS ))
  for ((g = G0; g < NG; g++)); do
    a=$((g*STEPS)); b=$((a + STEPS - 1)); [ $((g + 1)) -lt "$NG" ] || b=$KMAX
    echo "[$TAG] VTU group $g: steps $a-$b"
    # <tag>_s<k>.vtu and <tag>_t<k>.vtu by their step k (k = 100 t at dt 0.01); the relaxed frame counts as step 0; stall/tstop frames are left out
    ls "$D"/*.vtu 2>/dev/null | awk -v a=$a -v b=$b '{ f = $0
      if (f ~ /_relaxed\.vtu$/) k = 0
      else if (match(f, /_[st][0-9]+\.vtu$/)) k = substr(f, RSTART + 2, RLENGTH - 6) + 0
      else next
      if (k >= a && k <= b) print f }' | sort > "$LIST"
    part "$TAG: VTU steps $a-$b" || { echo "group $g failed: rerun with  bash push_run_data.sh $TAG vtu $g"; exit 1; }
  done
fi
echo "done: $REMOTE/$BR"
