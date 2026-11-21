#!/usr/bin/env python3
"""A1 comparison: the same 2D set-up on meshes of different element size.
For each case folder (remesh/a1/h<h>/ with out/ = run, grains/ = grain_post output, b1/ = b1_curvature output), over time:
  boundary length per misorientation bin (symmetry-reduced; <2, 2-5, 5-15, >=15 deg)  from grains/boundary_table_all.csv
  rotation front (um from the loaded edge)                                              from b1/b1_frames.csv
  E_stored, E_elastic, W_plate (work of the plate)                                      from out/*_log.csv
Prints, per quantity, the relative difference of each case to the finest one at a few times, and writes a1_compare.csv / a1_compare.png.
Pass criteria of the test plan: length per bin within 10 %, front within 0.1 um, plate work within 2 %.
usage: python3 remesh/a1_compare.py remesh/a1/h0.2 remesh/a1/h0.14 ... [--extra LABEL:DIR ...] [--times 1 2 3 5]
       (DIR of --extra: any run with out/ grains/ b1/ in the same layout, e.g. a production run's out_ + deliver_/grains)"""
import sys, os, glob, csv, argparse, re
import numpy as np
ap = argparse.ArgumentParser(); ap.add_argument("cases", nargs="+"); ap.add_argument("--times", type=float, nargs="+", default=[1, 2, 3, 4, 5])
ap.add_argument("--out", default="remesh/a1_compare")
a = ap.parse_args()
BINS = [(0, 2), (2, 5), (5, 15), (15, 1e9)]; BN = ["<2", "2-5", "5-15", ">=15"]
def lab(d):
    m = re.search(r"h([\d.]+)", os.path.basename(d.rstrip("/"))); return float(m.group(1)) if m else d
data = {}
for d in a.cases:
    h = lab(d); r = {}
    bt = os.path.join(d, "grains", "boundary_table_all.csv")
    if os.path.exists(bt):
        L = {}
        for row in csv.DictReader(open(bt)):
            t = float(row["t_ns"]); m = float(row["misorient_sym_deg"]); l = float(row["length_um"])
            v = L.setdefault(t, [0.0]*4)
            for i, (lo, hi) in enumerate(BINS):
                if lo <= m < hi: v[i] += l
        r["len"] = L
    b1 = os.path.join(d, "b1", "b1_frames.csv")
    if os.path.exists(b1):
        r["front"] = {float(x["t_ns"]): float(x["front_um"]) for x in csv.DictReader(open(b1))}
    lg = glob.glob(os.path.join(d, "out", "*_log.csv"))
    if lg:
        rows = list(csv.DictReader(open(lg[0])))
        r["log"] = {float(x["t_ns"]): x for x in rows}
    data[h] = r
hs = sorted([h for h in data if isinstance(h, float)])
ref = hs[0]
def at(series, t):
    if not series: return None
    ks = np.array(sorted(series)); k = ks[np.argmin(np.abs(ks - t))]
    return series[k] if abs(k - t) < 0.06 else None
out = []
print("reference: h = %g um (finest)" % ref)
for t in a.times:
    for h in hs:
        r, R = data[h], data[ref]; row = dict(t_ns=t, h_um=h)
        Lh, LR = at(r.get("len"), t), at(R.get("len"), t)
        if Lh and LR:
            for i, n in enumerate(BN): row["len_" + n] = Lh[i]; row["dlen_" + n] = (Lh[i] - LR[i])/max(LR[i], 1e-9)
        fh, fR = at(r.get("front"), t), at(R.get("front"), t)
        if fh is not None and fR is not None: row["front_um"] = fh; row["dfront_um"] = fh - fR
        gh, gR = at(r.get("log"), t), at(R.get("log"), t)
        if gh and gR:
            for c in ("E_stored", "E_elastic", "W_plate"):
                try:
                    x, y = float(gh[c]), float(gR[c]); row[c] = x; row["d" + c] = (x - y)/max(abs(y), 1e-12)
                except (KeyError, ValueError): pass
        out.append(row)
keys = sorted({k for r in out for k in r}, key=lambda k: (k not in ("t_ns", "h_um"), k))
os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
with open(a.out + ".csv", "w") as f:
    w = csv.DictWriter(f, fieldnames=keys); w.writeheader(); w.writerows(out)
for r in out:
    s = "t %4.1f h %5.3f |" % (r["t_ns"], r["h_um"])
    s += " ".join(" L%s %5.1f (%+5.0f%%)" % (n, r.get("len_" + n, np.nan), 100*r.get("dlen_" + n, np.nan)) for n in BN)
    s += " | front %5.2f (%+5.2f) | W_plate %+5.1f%%" % (r.get("front_um", np.nan), r.get("dfront_um", np.nan), 100*r.get("dW_plate", np.nan))
    print(s)
try:
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 5, figsize=(22, 4))
    for h in hs:
        L = data[h].get("len") or {}
        tt = sorted(L)
        for i in range(4): ax[i].plot(tt, [L[t][i] for t in tt], label="h %g" % h)
        F = data[h].get("front") or {}; tt = sorted(F); ax[4].plot(tt, [F[t] for t in tt], label="h %g" % h)
    for i in range(4): ax[i].set(title="boundary length, %s deg" % BN[i], xlabel="t (ns)", ylabel="um")
    ax[4].set(title="rotation front", xlabel="t (ns)", ylabel="um from loaded edge"); ax[0].legend()
    fig.tight_layout(); fig.savefig(a.out + ".png", dpi=110)
except Exception as e: print("plot failed:", e)
print("wrote", a.out + ".csv", a.out + ".png")
