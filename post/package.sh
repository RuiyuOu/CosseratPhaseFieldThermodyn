#!/bin/bash
# usage: bash package.sh <animout> <dest> [<datadir> <zipname>]   (paths may be relative; they are made absolute first)
set -e
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC="$(cd "$1" && pwd)"; mkdir -p "$2"; DST="$(cd "$2" && pwd)"
cd "$SRC"
for k in $(ls -d */ | tr -d /); do
  zip -qr "$DST/${k}_pngs.zip" "$k/png"
  python3 - "$k" "$DST" <<'P'
import sys,os,glob
from PIL import Image
k,dst=sys.argv[1:]; fs=sorted(glob.glob(f'{k}/png/*.png')); n=len(fs); pick=[fs[int(i*(n-1)/11)] for i in range(12)] if n>=12 else fs
ims=[Image.open(f) for f in pick]; w,h=ims[0].size; c=4; r=(len(ims)+c-1)//c
S=Image.new('RGB',(w*c,h*r),'white')
for i,im in enumerate(ims): S.paste(im,((i%c)*w,(i//c)*h))
S.save(f'{dst}/{k}_contact_sheet.png')
P
  cp "$k.gif" "$DST/"
done
cp "$HERE/SMOOTHING.md" "$DST/"
# raw data zip
if [ -n "$3" ]; then (cd "$3" && zip -qr "$DST/$4" . -x '*.out'); fi
