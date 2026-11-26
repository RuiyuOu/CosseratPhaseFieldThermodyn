#!/usr/bin/env python3
# kmt_patch.py -- run in the repo root.  Makes the Kocks-Mecking storage constant K and the annihilation distance d run-time options
# (-kmtk, -kmtd; defaults 0.05 and 5e-4 um = the old hard-coded values, so default runs are unchanged) in every cosserat_kwc*.edp, and
# forwards them from tag.sh as KMTK / KMTD only when they differ from the defaults (default PHYS keys stay byte-identical).
import glob, re, sys
OLD = 'real Kkmt = 0.05, dkmt = 5e-4;'
NEW = 'real Kkmt = atof(arg("-kmtk", "0.05")), dkmt = atof(arg("-kmtd", "5e-4"));'
done = 0
for f in sorted(glob.glob('cosserat_kwc*.edp')):
    s = open(f).read()
    if NEW in s: print(f, ': already patched'); done += 1; continue
    if s.count(OLD) != 1: print(f, ': WARNING no unique KMT line found, left untouched'); continue
    i = s.index(OLD); assert 'func string arg(' in s[:i], f + ': arg() must be defined before the KMT line'
    open(f, 'w').write(s.replace(OLD, NEW)); print(f, ': patched'); done += 1
assert done >= 1, 'no solver file patched'
t = open('tag.sh').read()
BLOCK = ('\n# Kocks-Mecking constants (solver defaults K 0.05, d 5e-4 um): appended only when changed, so the default key does not move\n'
         'if [ "${KMTK:-0.05}" != 0.05 ] || [ "${KMTD:-5e-4}" != 5e-4 ]; then PHYS="$PHYS -kmtk ${KMTK:-0.05} -kmtd ${KMTD:-5e-4}"; fi\n')
if '-kmtk' in t: print('tag.sh : already patched')
else:
    assert t.count('PHYS="$PHYS -diag ${DIAG:-1} -ckroll ${CKROLL:-20}"') == 1, 'tag.sh: last PHYS line not found'
    assert t.endswith('\n')
    open('tag.sh', 'w').write(t + BLOCK); print('tag.sh : patched')
print('done; check with: git diff --stat')
