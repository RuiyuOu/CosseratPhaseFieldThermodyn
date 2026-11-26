#!/usr/bin/env python3
"""panels4.py -- 2x2 panels (four times) with one shared colour bar on the right, for every field of the animations
(animate.py: phi, pressure, shear, ..., gnd) and the three grain figures (grain_theta, grain_territory, grain_misorientation).

usage (from post/):  BASE=<tag> python3 panels4.py ../out_<tag> ../anim_<tag> ../deliver_<tag>/panels4
env:  TIMES   "0 2.5 15 30"  target times in ns; the nearest output frame is used and its real time is printed in the panel
      FIELDS  optional subset, e.g. "phi gnd grain_misorientation"
      GRAD MERGE MINAREA SYM  grain segmentation as grain_post.py (defaults 2, 2, 0.03, 180)
      DPI     200
      GBOVER  "gnd rho_tot"  animation fields that get the grain boundaries (from the grain segmentation) drawn on top; "" for none
The animation colour limits are read from <anim_dir>/limits.json (the same scales as the GIFs); if it is missing they are computed as in
animate.py.  Grain theta colours use the limits of <deliver>/grains*/state.json if found, otherwise first-frame theta +- 20 deg; the
misorientation colour scale is 0 .. the largest misorientation in the four frames.  Output: <out_dir>/<field>.png and .pdf
"""
import os, sys, json, glob
from argparse import Namespace
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt, matplotlib.tri as mtri
from matplotlib.collections import LineCollection, PolyCollection
from matplotlib.colors import LogNorm, Normalize
from matplotlib.cm import ScalarMappable

D, ANIMDIR, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
sys.argv = ["animate.py", D, ANIMDIR]              # animate.py reads its run folder from argv at import
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import animate as A
import grain_count as gc

E = os.environ.get
TIMES = [float(x) for x in E("TIMES", "0 2.5 15 30").split()]
DPI = int(E("DPI", 200))
SEL = E("FIELDS", "").split()
PAR = Namespace(grad=float(E("GRAD", 2)), merge=float(E("MERGE", 2)), minarea=float(E("MINAREA", 0.03)), sym=float(E("SYM", 180)),
                lo=2.0, mid=5.0, hi=15.0, scan=False)
os.makedirs(OUT, exist_ok=True)

tt = np.array([t for t, _ in A.FR])
IDX = [int(np.argmin(np.abs(tt - x))) for x in TIMES]
print("frames used:", ", ".join("%g ns -> %.2f ns (%s)" % (x, A.FR[i][0], os.path.basename(A.FR[i][1])) for x, i in zip(TIMES, IDX)), flush=True)

# ---- colour limits of the animations
LF = os.path.join(ANIMDIR, "limits.json")
if os.path.exists(LF):
    A.LIM = {k: tuple(v) for k, v in json.load(open(LF)).items()}
    new = [k for k in A.ANIMS if k not in A.LIM]
    if new: A.LIM.update(A.limits(new))
    A.MAG = A.LIM.get("_mag", (A.MAG,))[0]
else:
    print("no %s: computing colour limits as animate.py does" % LF, flush=True)
    A.LIM = A.limits(); A.LIM["_mag"] = (A.MAG, A.MAG)


def grid():
    fig, ax = plt.subplots(2, 2, figsize=(10.2, 9.0), dpi=DPI, constrained_layout=True)
    return fig, ax.ravel()


def stamp(x, t):
    x.text(0.02, 0.97, "t = %.1f ns" % t, transform=x.transAxes, va="top", fontsize=10, color="w", bbox=dict(fc="k", alpha=.55, lw=0))


def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(OUT, "%s.%s" % (name, ext)))
    plt.close(fig); print("written", name, flush=True)


# ---- grain segmentation of the four frames (used by the grain figures and the boundary overlay)
_G = {}
def segs(i):
    if i not in _G:
        r = gc.analyse_frame(A.FR[i][1], PAR)
        if r is None: sys.exit("grain segmentation failed for %s" % A.FR[i][1])
        _G[i] = r
    return _G[i]
GBOVER = E("GBOVER", "gnd rho_tot").split()


def gb_overlay(x, i):
    M, seg, rows, bnd = segs(i)
    L = [M.P[M.en[e]] for b in bnd for e in b["edges"]]
    x.add_collection(LineCollection(L, colors="k", linewidths=1.8, alpha=0.7))
    x.add_collection(LineCollection(L, colors="#00e5ff", linewidths=0.9))


# ---- animation fields (plus the total dislocation density rho_S + rho_G, m^-2)
O = {i: A.get(i) for i in IDX}
A.ANIMS["rho_tot"] = ("Total dislocation density $\\rho_S+\\rho_G$", "rho_tot", "magma", False, True, "m$^{-2}$")
for i in IDX:
    d = A.w.load(A.FR[i][1], A.M); g, s_ = np.asarray(d["rhoG"], float), np.asarray(d["rhoS"], float)
    if len(g) != len(s_):
        g, s_ = [x[A.M.tri].mean(1) if len(x) == len(A.M.p) else x for x in (g, s_)]
    O[i][1]["rho_tot"] = (g + s_) * 1e12
if "rho_tot" not in A.LIM:
    LG = os.path.join(ANIMDIR, "limits_gb.json")
    if os.path.exists(LG) and "rhotot_gb" in json.load(open(LG)): A.LIM["rho_tot"] = tuple(json.load(open(LG))["rhotot_gb"])
    else:
        a = np.concatenate([O[i][1]["rho_tot"] for i in IDX]); a = a[a > 0]
        A.LIM["rho_tot"] = (max(float(np.percentile(a, 2)), 1e8), max(float(np.percentile(a, 99.8)), 1e9))
for k, (ttl, key, cm, sym, lg, un) in A.ANIMS.items():
    if SEL and k not in SEL: continue
    lo, hi = A.LIM[k]
    fig, axs = grid(); pc = None
    for x, i in zip(axs, IDX):
        t, o = O[i]; v = o[key]
        P = A.M.p.copy()
        if k == "deformation": P = P + A.MAG * o["u"]
        tr = mtri.Triangulation(P[:, 0], P[:, 1], A.M.tri)
        if k == "lattice_slip": v = ((v + 90.0) % 180.0) - 90.0
        if lg: v = np.maximum(v, lo)
        kw = dict(cmap=cm, norm=LogNorm(lo, hi)) if lg else dict(cmap=cm, vmin=lo, vmax=hi)
        pc = x.tripcolor(tr, v, **kw, shading="flat") if len(v) == A.M.nt else x.tripcolor(tr, v, **kw, shading="gouraud")
        if k == "lattice_slip": A.slip_overlay(x, o)
        if k in GBOVER: gb_overlay(x, i)
        if k == "deformation":
            x.plot([0, 10, 10, 0, 0], [0, 0, 10, 10, 0], "w--", lw=.8); x.set_xlim(-1.2, 11.2); x.set_ylim(-1.2, 11.2)
        else:
            x.set_xlim(A.M.p[:, 0].min(), A.M.p[:, 0].max()); x.set_ylim(A.M.p[:, 1].min(), A.M.p[:, 1].max())
        x.set_aspect("equal"); x.set_xlabel("x (µm)"); x.set_ylabel("y (µm)"); stamp(x, t)
    if k == "deformation": ttl = "Deformation $|u|$, displacement magnified %.0fx (dashed = undeformed)" % A.MAG
    if k == "gnd": ttl = "GND density (Cermelli–Gurtin): $|G|/b$, $G=F^p\\,\\mathrm{Curl}\\,F^p$"
    if k in GBOVER: ttl += "\nlines: boundaries between grain territories (MINAREA %g, merge %g$^\\circ$)" % (PAR.minarea, PAR.merge)
    fig.suptitle(ttl, fontsize=11)
    fig.colorbar(pc, ax=list(axs), shrink=0.9, label=un)
    save(fig, k)

# ---- grain fields
GK = ["grain_theta", "grain_territory", "grain_misorientation"]
if not SEL or any(k in SEL for k in GK):
    lim = None
    for sf in sorted(glob.glob(os.path.join(os.path.dirname(OUT.rstrip("/")), "grains*", "state.json"))):
        st = json.load(open(sf))
        if st.get("lim"): lim = st["lim"]; print("grain theta limits from", sf, flush=True); break
    if lim is None:
        _, _, th0, _, _ = gc.load(A.FR[0][1]); lim = [float(th0.min()) - 20.0, float(th0.max()) + 20.0, 150.0]
    th0, th1 = lim[0], lim[1]
    G = {i: segs(i) for i in IDX}
    mmax = max([b["misorient_deg"] for i in IDX for b in G[i][3]] + [1.0])

    def frame_ax(x, M, t):
        x.set_aspect("equal"); x.set_xlim(M.P[:, 0].min(), M.P[:, 0].max()); x.set_ylim(M.P[:, 1].min(), M.P[:, 1].max())
        x.set_xlabel("x (µm)"); x.set_ylabel("y (µm)"); stamp(x, t)

    if not SEL or "grain_theta" in SEL:
        fig, axs = grid()
        for x, i in zip(axs, IDX):
            M = G[i][0]; tr = mtri.Triangulation(M.P[:, 0], M.P[:, 1], M.T)
            pc = x.tripcolor(tr, M.th, shading="gouraud", cmap="twilight", vmin=th0, vmax=th1); frame_ax(x, M, A.FR[i][0])
        fig.suptitle("Lattice orientation $\\theta_{lat}$", fontsize=11); fig.colorbar(pc, ax=list(axs), shrink=0.9, label="deg"); save(fig, "grain_theta")

    if not SEL or "grain_territory" in SEL:
        fig, axs = grid(); cm = plt.get_cmap("twilight")
        for x, i in zip(axs, IDX):
            M, seg, rows, bnd = G[i]
            mug = seg["mu"][np.maximum(seg["terr"], 0)]; col = cm(np.clip((mug - th0) / (th1 - th0), 0, 1))
            col[seg["plat"] < 0] = col[seg["plat"] < 0] * 0.35 + 0.65
            x.add_collection(PolyCollection(M.P[M.T], facecolors=col, edgecolors=col, linewidths=0.3))
            for g in rows:
                x.text(g["cx"], g["cy"], str(g["id"]), ha="center", va="center", fontsize=7, weight="bold", bbox=dict(fc="w", ec="none", alpha=0.55, pad=0.5))
            frame_ax(x, M, A.FR[i][0]); x.set_title("%d grains" % seg["G"], fontsize=9)
        fig.suptitle("Grain territories: plateaus (solid), ramps and cores (pale); colour = grain mean $\\theta$", fontsize=11)
        fig.colorbar(ScalarMappable(norm=Normalize(th0, th1), cmap=cm), ax=list(axs), shrink=0.9, label="grain mean $\\theta$ (deg)")
        save(fig, "grain_territory")

    if not SEL or "grain_misorientation" in SEL:
        fig, axs = grid(); lc = None
        for x, i in zip(axs, IDX):
            M, seg, rows, bnd = G[i]
            x.add_collection(PolyCollection(M.P[M.T], facecolors=plt.get_cmap("Greys")(0.08 + 0.12 * (seg["plat"] >= 0)), edgecolors="none"))
            segs, cc = [], []
            for b in bnd:
                for e in b["edges"]:
                    segs.append(M.P[M.en[e]]); cc.append(b["misorient_deg"])
            lc = LineCollection(segs, array=np.array(cc), cmap="turbo", linewidths=1.6); lc.set_clim(0, mmax); x.add_collection(lc)
            for b in bnd:
                m = M.emid[b["edges"]].mean(0)
                x.text(m[0], m[1], "%.1f" % b["misorient_deg"], fontsize=5.5, ha="center", va="center", bbox=dict(fc="w", ec="none", alpha=0.7, pad=0.3))
            frame_ax(x, M, A.FR[i][0])
        fig.suptitle("Boundaries between grain territories (label = misorientation, deg)", fontsize=11)
        fig.colorbar(lc, ax=list(axs), shrink=0.9, label="misorientation (deg, raw)"); save(fig, "grain_misorientation")
print("done:", OUT)
