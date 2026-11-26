#!/bin/bash
# push_chunked.sh -- push whole folders (VTUs included) to the data branch in several pushes of at most CHUNKMB each,
# one after another, so that no single push exceeds GitHub's push size limit.
#   bash push_chunked.sh "<commit message>" <path> [<path> ...]
#   env: CHUNKMB (1500)  MAXMB (95: larger files are skipped and listed)  DRY=1 (only report sizes)
# Each chunk is committed and pushed by push_paths.sh; a rerun skips nothing but re-pushes nothing new (safe to resume).
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
MSG=${1:?usage: bash push_chunked.sh "<message>" <path> [...]}; shift
CH=$(( ${CHUNKMB:-1500} * 1024 * 1024 )); MX=$(( ${MAXMB:-95} * 1024 * 1024 ))
LIST=$(mktemp -d)/list; trap 'rm -rf "$(dirname "$LIST")"' EXIT
for p in "$@"; do
  [ -e "$p" ] || { echo "  missing: $p" >&2; continue; }
  find "$p" -type f ! -name '*.tmp' ! -name '*.lock' -printf '%s\t%p\n'
done | sort -t$'\t' -k2 > "$LIST"
awk -F'\t' -v mx=$MX '$1>mx{print "  SKIPPED (too large): " $2}' "$LIST" >&2
awk -F'\t' -v mx=$MX -v ch=$CH -v out="$LIST" 'BEGIN{k=1;s=0} $1<=mx{if(s+$1>ch && s>0){k++;s=0} s+=$1; print $2 > (out "." k)}' "$LIST"
N=$(ls "$LIST".* 2>/dev/null | wc -l)
echo "$(awk -F'\t' -v mx=$MX '$1<=mx{n++;s+=$1}END{printf "%d files, %.2f GB", n, s/2^30}' "$LIST") in $N pushes"
[ "${DRY:-0}" = 1 ] && exit 0
for k in $(seq 1 "$N"); do
  echo "--- push $k/$N"; mapfile -t F < "$LIST.$k"
  VTU=1 MAXMB=${MAXMB:-95} bash push_paths.sh "$MSG ($k/$N)" "${F[@]}"
done
