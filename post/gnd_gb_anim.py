#!/usr/bin/env python3
"""gnd_gb_anim.py -- animations of dislocation densities with the boundaries between grain territories drawn on top in red:
  gnd_gb      GND density rho_G (the solver's Cermelli-Gurtin field, as the GND animation)
  rhotot_gb   total dislocation density rho_S + rho_G

usage (from post/):  BASE=<tag> NP=4 python3 gnd_gb_anim.py ../out_<tag> ../anim_<tag> [../deliver_<tag>]
writes, for both names, <deliver_dir>/<name>/png/<name>_NNNN.png (one per output frame), <deliver_dir>/<name>/<name>.gif and
<deliver_dir>/<name>/<name>_contact_sheet.png (12 evenly spaced frames, 3 x 4), and a copy of the GIF to <anim_dir>/<name>.gif.
<deliver_dir> defaults to deliver_<tag> next to the run folder. Frames drawn by an earlier version in <anim_dir>/<name>/png are moved over.
env:  WHICH  "gnd_gb rhotot_gb" (default both)
      GRAD MERGE MINAREA SYM   grain segmentation as grain_post.py (defaults 2, 2, 0.03, 180)
      EVERY  use every n-th output frame (default 1);  NP worker processes (default 4);  FORCE=1 redraws existing frames;  DPI (110)
Colour scales (log): GND from <anim_dir>/limits.json (the GND animation's scale) if present; total density from 12 evenly spaced frames
(2nd and 99.8th percentiles), stored in <anim_dir>/limits_gb.json and reused on later calls so a resumed series keeps its scale.
"""
import os, sys, json, shutil
from argparse import Namespace
from multiprocessing import Pool
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt, matplotlib.tri as mtri
from matplotlib.collections import LineCollection
from matplotlib.colors import LogNorm

D, ANIMDIR = sys.argv[1], sys.argv[2]
_b = os.path.basename(os.path.normpath(D))
DELDIR = sys.argv[3] if len(sys.argv) > 3 else os.path.join(os.path.dirname(os.path.normpath(D)), "deliver_" + (_b[4:] if _b.startswith("out_") else _b))
sys.argv = ["animate.py", D, ANIMDIR]              # animate.py reads its run folder from argv at import
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import animate as A
import grain_count as gc
from PIL import Image

E = os.environ.get
PAR = Namespace(grad=float(E("GRAD", 2)), merge=float(E("MERGE", 2)), minarea=float(E("MINAREA", 0.03)), sym=float(E("SYM", 180)),
                lo=2.0, mid=5.0, hi=15.0, scan=False)
EVERY, DPI, FORCE = int(E("EVERY", 1)), int(E("DPI", 110)), E("FORCE", "0") == "1"
GBC = {"gnd_gb": "red", "rhotot_gb": "lime"}; GBN = {"gnd_gb": "red", "rhotot_gb": "green"}
WHICH = E("WHICH", "gnd_gb rhotot_gb").split()
SPEC = {"gnd_gb": ("GND density $\\rho_G=|G|/b$", "magma"),
        "rhotot_gb": ("total dislocation density $\\rho_S+\\rho_G$", "magma")}


def fields(i):
    """GND and total density (m^-2) of frame i, element-wise when both fields are element fields"""
    d = A.w.load(A.FR[i][1], A.M)
    g, s = np.asarray(d["rhoG"], float), np.asarray(d["rhoS"], float)
    nt, nn = A.M.nt, len(A.M.p)
    to_cell = lambda x: x[A.M.tri].mean(1) if len(x) == nn else x
    if len(g) != len(s): g, s = to_cell(g), to_cell(s)
    return {"gnd_gb": g * 1e12, "rhotot_gb": (g + s) * 1e12}


# ---- colour limits
LG = os.path.join(ANIMDIR, "limits_gb.json"); LIMS = json.load(open(LG)) if os.path.exists(LG) else {}
LF = os.path.join(ANIMDIR, "limits.json")
if "gnd_gb" not in LIMS:
    if os.path.exists(LF) and "gnd" in json.load(open(LF)): LIMS["gnd_gb"] = list(json.load(open(LF))["gnd"])
if any(k not in LIMS for k in WHICH):
    acc = {k: [] for k in SPEC}
    for i in np.unique(np.linspace(0, len(A.FR) - 1, min(len(A.FR), 12)).astype(int)):
        f = fields(i)
        for k in SPEC: acc[k].append(f[k])
    for k in SPEC:
        if k in LIMS: continue
        a = np.concatenate(acc[k]); a = a[a > 0]
        LIMS[k] = [max(float(np.percentile(a, 2)), 1e8), max(float(np.percentile(a, 99.8)), 1e9)] if len(a) else [1e8, 1e16]
os.makedirs(ANIMDIR, exist_ok=True); json.dump(LIMS, open(LG, "w"), indent=1)


def png(name, i):
    return os.path.join(DELDIR, name, "png", "%s_%04d.png" % (name, i))


def old_png(name, i):
    return os.path.join(ANIMDIR, name, "png", "%s_%04d.png" % (name, i))


def contact_sheet(name, fs, out, n=12, cols=4):
    pick = [fs[int(round(j))] for j in np.linspace(0, len(fs) - 1, min(n, len(fs)))]
    ims = [Image.open(f).convert("RGB") for f in pick]
    w, h = ims[0].size; rows = (len(ims) + cols - 1) // cols
    sheet = Image.new("RGB", (cols * w, rows * h), "white")
    for j, im in enumerate(ims): sheet.paste(im.resize((w, h)), ((j % cols) * w, (j // cols) * h))
    sheet.save(out)


def render(i):
    try:
        need = [k for k in WHICH if FORCE or not os.path.exists(png(k, i))]
        if not need: return i
        t = A.FR[i][0]; F = fields(i)
        r = gc.analyse_frame(A.FR[i][1], PAR); n = 0; L = []
        if r is not None:
            M, seg, rows, bnd = r; n = seg["G"]; L = [M.P[M.en[e]] for b in bnd for e in b["edges"]]
        tr = mtri.Triangulation(A.M.p[:, 0], A.M.p[:, 1], A.M.tri)
        for k in need:
            lo, hi = LIMS[k]; v = np.maximum(F[k], lo)
            fig, ax = plt.subplots(figsize=(5.6, 5.0), dpi=DPI)
            pc = ax.tripcolor(tr, v, cmap=SPEC[k][1], norm=LogNorm(lo, hi), shading="flat" if len(v) == A.M.nt else "gouraud")
            if L: ax.add_collection(LineCollection(L, colors=GBC[k], linewidths=1.2))
            ax.set_aspect("equal"); ax.set_xlim(A.M.p[:, 0].min(), A.M.p[:, 0].max()); ax.set_ylim(A.M.p[:, 1].min(), A.M.p[:, 1].max())
            ax.set_xlabel("x (µm)"); ax.set_ylabel("y (µm)")
            ax.set_title("%s with grain boundaries (%s, %d grains)" % (SPEC[k][0], GBN[k], n), fontsize=9)
            ax.text(0.02, 0.97, "t = %.2f ns" % t, transform=ax.transAxes, va="top", fontsize=9, color="w", bbox=dict(fc="k", alpha=.5, lw=0))
            cb = fig.colorbar(pc, ax=ax, shrink=.85); cb.set_label("m$^{-2}$")
            fig.tight_layout(); fig.savefig(png(k, i)); plt.close(fig)
        return i
    except Exception as e:                        # one bad frame must not stop the series
        print("gnd_gb_anim: frame %d failed: %s" % (i, e), flush=True)
        return None


if __name__ == "__main__":
    for k in WHICH: os.makedirs(os.path.dirname(png(k, 0)), exist_ok=True)
    idx = list(range(0, len(A.FR), EVERY))
    if idx[-1] != len(A.FR) - 1: idx.append(len(A.FR) - 1)
    for k in WHICH:                               # frames from the earlier version (written under the anim folder)
        nm = 0
        for i in idx:
            if os.path.exists(old_png(k, i)) and not os.path.exists(png(k, i)): shutil.move(old_png(k, i), png(k, i)); nm += 1
        if nm: print("moved %d earlier %s frames to %s" % (nm, k, os.path.dirname(png(k, 0))), flush=True)
    todo = [i for i in idx if any(FORCE or not os.path.exists(png(k, i)) for k in WHICH)]
    print("gnd_gb_anim: %s, %d frames, %d to render (MINAREA %g, merge %g, SYM %g)" % (" ".join(WHICH), len(idx), len(todo), PAR.minarea, PAR.merge, PAR.sym), flush=True)
    if todo:
        with Pool(int(E("NP", 4))) as p:
            p.map(render, todo, chunksize=1)
    for k in WHICH:
        fs = [png(k, i) for i in idx if os.path.exists(png(k, i))]
        if fs:
            ims = [Image.open(f).convert("P", palette=Image.ADAPTIVE) for f in fs]
            gif = os.path.join(DELDIR, k, k + ".gif")
            ims[0].save(gif, save_all=True, append_images=ims[1:], duration=80, loop=0)
            shutil.copy2(gif, os.path.join(ANIMDIR, k + ".gif"))
            cs = os.path.join(DELDIR, k, k + "_contact_sheet.png"); contact_sheet(k, fs, cs)
            print("written", gif, "and", cs, "(%d frames)" % len(fs), flush=True)
        else:
            print("no %s frames written: see the 'frame N failed' lines above" % k, flush=True)
