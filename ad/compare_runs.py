#!/usr/bin/env python3
"""compare_runs.py RUN_A RUN_B  -- per-step comparison of two runs of the same set-up (logs, last common VTU, wall times)."""
import sys, os, csv, glob, re
A, B = sys.argv[1], sys.argv[2]
def log(d):
    f = glob.glob(os.path.join(d, '*_log.csv'))[0]; return {r['step']: r for r in csv.DictReader(open(f))}
a, b = log(A), log(B)
steps = [s for s in a if s in b]
print("%d common logged steps" % len(steps))
cols = ['iters', 'E_elastic', 'KE_trans', 'Q_heat', 'max_abs_slip', 'max_dT', 'max_abs_dtheta_deg', 'max_rhoG_um2', 'min_phi', 'audit_full']
worst = {c: 0. for c in cols}; itd = 0; ita = itb = 0
for s in steps:
    for c in cols:
        x, y = a[s].get(c, 'nan'), b[s].get(c, 'nan')
        if 'nan' in (x, y): continue
        x, y = float(x), float(y)
        if c == 'iters': ita += x; itb += y; itd += (x != y); continue
        worst[c] = max(worst[c], abs(x - y)/max(abs(x), 1e-12))
print("Newton iterations: total %d (A) vs %d (B), %d steps differ" % (ita, itb, itd))
for c in cols[1:]: print("  max relative difference %-20s %.2e" % (c, worst[c]))
for d in (A, B):
    of = glob.glob(os.path.join(d, '*.out'))
    w = [l for l in open(of[0]) if 'wall time' in l] if of else []
    print("%s: %s" % (os.path.basename(d), (w[-1].strip()[:110] if w else 'no wall time line')))
try:
    import meshio, numpy as np
    va = {os.path.basename(f).split('_s')[-1]: f for f in glob.glob(os.path.join(A, '*_s*.vtu'))}
    vb = {os.path.basename(f).split('_s')[-1]: f for f in glob.glob(os.path.join(B, '*_s*.vtu'))}
    com = sorted(set(va) & set(vb), key=lambda k: int(re.sub(r'\D', '', k) or 0))
    if com:
        k = com[-1]; ma, mb = meshio.read(va[k]), meshio.read(vb[k])
        print("last common VTU s%s:" % k.replace('.vtu', ''))
        for f in ['u', 'phi', 'T']:
            x, y = ma.point_data[f], mb.point_data[f]; print("  %-8s max |A-B| %.2e  (max |A| %.2e)" % (f, abs(x - y).max(), abs(x).max()))
        for f in ['p1', 'sig_11', 'rhoS', 'rhoG']:
            x, y = ma.cell_data[f][0], mb.cell_data[f][0]; print("  %-8s max |A-B| %.2e  (max |A| %.2e)" % (f, abs(x - y).max(), abs(x).max()))
except ImportError:
    print("(meshio not available: VTU comparison skipped; source ~/dynpost_venv/bin/activate)")
