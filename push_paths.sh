#!/bin/bash
# push_paths.sh -- commit arbitrary files or folders to the data branch and push them, without touching the checked-out branch.
#   bash push_paths.sh "<commit message>" <path> [<path> ...]
#   env: BRANCH (data)  REMOTE (origin)  MAXMB (95: larger files are skipped)
# Folders are taken recursively; VTUs are left out unless VTU=1. Same separate-index method as push_run_data.sh.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
MSG=${1:?usage: bash push_paths.sh "<message>" <path> [...]}; shift
[ $# -gt 0 ] || { echo "no paths given"; exit 1; }
BR=${BRANCH:-data}; REMOTE=${REMOTE:-origin}; MAXMB=${MAXMB:-95}
git rev-parse --git-dir >/dev/null 2>&1 || { echo "not inside the git repository"; exit 1; }
CUR=$(git symbolic-ref --short -q HEAD || true)
[ "$BR" != "$CUR" ] || { echo "BRANCH=$BR is checked out here: use another branch"; exit 1; }
export GIT_INDEX_FILE="$(git rev-parse --absolute-git-dir)/index.paths.$$"
LIST=$(mktemp); trap 'rm -f "$GIT_INDEX_FILE" "$LIST"' EXIT
if git rev-parse -q --verify "refs/heads/$BR" >/dev/null; then P=$(git rev-parse "refs/heads/$BR")
elif git ls-remote --exit-code --heads "$REMOTE" "$BR" >/dev/null 2>&1; then git fetch -q "$REMOTE" "$BR:refs/heads/$BR"; P=$(git rev-parse "refs/heads/$BR")
else P=""; fi
if [ -n "$P" ]; then git read-tree "$P"; else git read-tree --empty; fi
for p in "$@"; do
  [ -e "$p" ] || { echo "  missing: $p" >&2; continue; }
  if [ "${VTU:-0}" = 1 ]; then find "$p" -type f ! -name '*.tmp' ! -name '*.lock'; else find "$p" -type f ! -name '*.vtu' ! -name '*.tmp' ! -name '*.lock'; fi
done | sort | while IFS= read -r f; do
  if [ "$(stat -c %s "$f")" -gt $((MAXMB*1024*1024)) ]; then echo "  SKIPPED (> $MAXMB MB): $f" >&2; else echo "$f"; fi
done > "$LIST"
n=$(wc -l < "$LIST"); [ "$n" -gt 0 ] || { echo "no files to push"; exit 0; }
echo "$n files, $(xargs -a "$LIST" -d '\n' du -ch | tail -1 | cut -f1)"
xargs -a "$LIST" -d '\n' git update-index --add --
T=$(git write-tree)
if [ -n "$P" ] && [ "$T" = "$(git rev-parse "$P^{tree}")" ]; then echo "already on $BR: nothing to commit"; else
  C=$(git commit-tree "$T" ${P:+-p "$P"} -m "$MSG"); git update-ref "refs/heads/$BR" "$C" ${P:+"$P"}; echo "commit ${C:0:12}"
fi
git -c http.postBuffer=524288000 push -q "$REMOTE" "refs/heads/$BR:refs/heads/$BR" && echo "pushed to $REMOTE/$BR"
