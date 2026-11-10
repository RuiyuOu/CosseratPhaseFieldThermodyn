#!/usr/bin/env python3
"""grain_post.py  --  grain count and boundary-misorientation history of a run, incremental.

usage:  BASE=<tag> NP=2 python3 grain_post.py <run_dir> <out_dir> [<anim_dir>]

Takes the same frames as animate.py (wavepost.frame_list: the relaxed frame + every frame of <BASE>.pvd) and runs
grain_count.py on each.  Only frames that are not yet in <out_dir>/grain_history.csv are analysed, so a second call after
the run has continued (or after a crash and resume) costs only the new frames.  A change of the grain parameters
re-analyses everything (FORCE=1 does too).

env:  GRAD (2 deg/um)  MERGE (2 deg)  MINAREA (0.1 um^2)  SYM (90 deg; 360 = none)  LO/MID/HI (2/5/15 deg bin edges)
      NP (processes)  FORCE=1
anim_dir (optional, normally ../anim_<tag>): per-frame figures and GIFs grain_theta, grain_territory, grain_misorientation are written there as
          <anim_dir>/<name>/png/<name>_<frame>.png and <anim_dir>/<name>.gif, so package.sh makes their contact sheets / PNG zips like the other
          animations.  Fixed colour scales (theta range of the first frame +-20 deg, misorientation 0..150 deg); territories are coloured by
          their mean theta.  Incremental: only frames without a figure are drawn, a GIF is rebuilt only when frames were added.
out_dir:  grain_table_all.csv     one row per grain and frame (t, step, file, id, areas, centre, theta, orig_label, neighbours)
          boundary_table_all.csv  one row per neighbour pair and frame (grain ids and thetas, misorientation raw / mod SYM, length, width, mean phi, midpoint)
          grain_history.csv       one row per frame: t, step, grains, boundary length, % of length per misorientation bin (raw, mod SYM)
          grain_history.png       number of grains and the bin fractions against time
          grains_first_* / grains_last_*   full per-grain / per-boundary tables, histogram and the three figures (theta, territory, misorientation) of the first (relaxed) and last frame
          state.json              parameters and the frames done (do not edit)
"""
import sys, os, re, json, csv, time
from argparse import Namespace
from multiprocessing import Pool
import numpy as np
import wavepost as w
import grain_count as gc

D, OUT = sys.argv[1], sys.argv[2]
ANIM = sys.argv[3] if len(sys.argv) > 3 else None
FIGS = ("theta", "territory", "misorientation")
LIM = None
BASE = os.environ.get("BASE", "A_dyn")
E = os.environ.get
PAR = Namespace(grad=float(E("GRAD", 2)), merge=float(E("MERGE", 2)), minarea=float(E("MINAREA", 0.1)), sym=float(E("SYM", 90)),
                lo=float(E("LO", 2)), mid=float(E("MID", 5)), hi=float(E("HI", 15)), scan=False)
KEYS = ["grad", "merge", "minarea", "sym", "lo", "mid", "hi"]


def step_of(fn):
    m = re.search(r"_s(\d+)\.vtu$", fn)
    return int(m.group(1)) if m else -1


def fig_paths(i):
    return {k: os.path.join(ANIM, "grain_" + k, "png", "grain_%s_%04d.png" % (k, i)) for k in FIGS}


def work(args):
    i, t, fn, need_row, need_figs = args
    try:
        r = gc.analyse_frame(fn, PAR)
        if r is None:
            return None
        row = tabs = None
        if need_row:
            row = gc.summary_from(*r, PAR); row.update(t_ns=t, step=step_of(fn), file=os.path.basename(fn))
            M_, seg_, rows_, bnd_ = r; nn = [0] * seg_["G"]
            for b in bnd_: nn[b["i"]] += 1; nn[b["j"]] += 1
            key = dict(t_ns=t, step=step_of(fn), file=os.path.basename(fn))
            gr = [dict(key, id=g["id"], plateau_area_um2=g["plateau_area"], territory_area_um2=g["territory_area"], cx_um=g["cx"], cy_um=g["cy"],
                       theta_deg=g["theta_deg"], theta_std_deg=g["theta_std_deg"], theta_sym_deg=g["theta_sym_deg"], orig_label=g["orig_label"],
                       n_neighbours=nn[g["id"]]) for g in rows_]
            br = []
            for b in bnd_:
                m_ = M_.emid[b["edges"]].mean(0)
                br.append(dict(key, grain_i=b["i"], grain_j=b["j"], theta_i_deg=seg_["mu"][b["i"]], theta_j_deg=seg_["mu"][b["j"]],
                               misorient_deg=b["misorient_deg"], misorient_sym_deg=b["misorient_sym_deg"], length_um=b["length_um"],
                               width_um=b["width_um"], mean_phi_at_edge=b["phi_at_edge"], mid_x_um=m_[0], mid_y_um=m_[1]))
            tabs = (gr, br)
        if need_figs and ANIM:
            M, seg, rows, bnd = r; P = fig_paths(i)
            for p in P.values(): os.makedirs(os.path.dirname(p), exist_ok=True)
            gc.plot(M, seg, rows, bnd, P, PAR, lim=LIM, t=t)
        return (row, tabs) if row is not None else None
    except Exception as e:                       # one bad frame must not stop the series
        print("grain_post: %s failed: %s" % (os.path.basename(fn), e), flush=True)
        return None


def init(l):
    global LIM
    LIM = l


def detail(fn, pfx):
    a = Namespace(**vars(PAR)); a.vtu = fn; a.out = pfx
    import io, contextlib
    with contextlib.redirect_stdout(io.StringIO()):
        gc.report(a, pfx)


def plot(rows, fn):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    t = np.array([r["t_ns"] for r in rows]); names = [n for n, _, _ in gc.bins_of(PAR)]
    fig, ax = plt.subplots(3, 1, figsize=(8, 10), sharex=True)
    ax[0].plot(t, [r["n_grains"] for r in rows], "o-", ms=3, label="all plateaus")
    ax[0].plot(t, [r["n_ge1um2"] for r in rows], "s-", ms=3, label="plateau area >= 1 um$^2$")
    ax[0].set_ylabel("number of grains"); ax[0].legend(); ax[0].grid(alpha=.3)
    ax[0].set_title("grain-count history (grad %g deg/um, merge %g deg, min area %g um$^2$)" % (PAR.grad, PAR.merge, PAR.minarea), fontsize=10)
    for k, tag in ((1, "raw"), (2, "sym")):
        y = np.array([[r[f"pct_{tag}_{n}"] for n in names] for r in rows]).T
        ax[k].stackplot(t, y, labels=names, colors=["#bbbbbb", "#4c9be8", "#f5a623", "#d0021b"])
        ax[k].set_ylim(0, 100); ax[k].set_ylabel("% of boundary length")
        ax[k].set_title("misorientation = " + ("raw difference" if tag == "raw" else "reduced mod %g deg" % PAR.sym), fontsize=9)
    ax[1].legend(loc="center left", bbox_to_anchor=(1.0, 0.5), title="deg", fontsize=8)
    ax[2].set_xlabel("t (ns)"); plt.tight_layout(); plt.savefig(fn, dpi=100); plt.close()


if __name__ == "__main__":
    t0 = time.time()
    os.makedirs(OUT, exist_ok=True)
    FR = w.frame_list(D, BASE)
    if not FR:
        sys.exit("grain_post: no frames in %s (base %s)" % (D, BASE))
    sf = os.path.join(OUT, "state.json"); hf = os.path.join(OUT, "grain_history.csv")
    par = {k: getattr(PAR, k) for k in KEYS}
    old = []; params_changed = False
    if os.path.exists(hf) and os.path.exists(sf) and E("FORCE", "0") != "1":
        st = json.load(open(sf))
        if st.get("params") == par:
            old = [dict(r) for r in csv.DictReader(open(hf))]
            for r in old:
                for k in r:
                    if k != "file": r[k] = float(r[k])
        else:
            print("grain_post: parameters changed, re-analysing all frames"); params_changed = True
    gtf = os.path.join(OUT, "grain_table_all.csv"); btf = os.path.join(OUT, "boundary_table_all.csv")
    def read_tab(fn):
        return [dict(r) for r in csv.DictReader(open(fn))] if (old and os.path.exists(fn)) else []
    gtab, btab = read_tab(gtf), read_tab(btf)
    tab_files = {r["file"] for r in gtab}                  # frames whose per-grain table exists (older output folders have none)
    have = {r["file"] for r in old}
    # colour scales of the animation frames: fixed from the first frame of the first call (kept in state.json)
    stp = json.load(open(sf)) if os.path.exists(sf) else {}
    lim = stp.get("lim") if (stp.get("params") == par and E("FORCE", "0") != "1") else None
    if lim is None:
        pts0, tri0, th0, _, _ = gc.load(FR[0][1])
        lim = [float(th0.min()) - 20.0, float(th0.max()) + 20.0, 150.0]
    # a frame needs work if it has no history row, or (with an anim dir) a missing figure; a parameter change makes `old` empty
    def figs_ok(i):
        return (not ANIM) or all(os.path.exists(p) for p in fig_paths(i).values())
    todo = []
    for i, (t, f) in enumerate(FR):
        nr = (os.path.basename(f) not in have) or (os.path.basename(f) not in tab_files); nf = not figs_ok(i)
        if params_changed or E("FORCE", "0") == "1": nf = bool(ANIM)
        if nr or nf: todo.append((i, t, f, nr, nf))
    print("grain_post: %d frames, %d already done, %d to analyse%s" % (len(FR), len(FR) - len(todo), len(todo), " (figures into %s)" % ANIM if ANIM else ""), flush=True)
    new = []
    if todo:
        with Pool(int(E("NP", 2)), initializer=init, initargs=(lim,)) as p:
            res = [r for r in p.imap(work, todo, chunksize=1) if r is not None]
        new = [r[0] for r in res]
        redone = {r[0]["file"] for r in res}
        gtab = [r for r in gtab if r["file"] not in redone] + [x for r in res for x in r[1][0]]
        btab = [r for r in btab if r["file"] not in redone] + [x for r in res for x in r[1][1]]
    if ANIM and todo and any(x[4] for x in todo):
        from PIL import Image
        for k in FIGS:
            fs = [fig_paths(i)[k] for i in range(len(FR))]
            fs = [f for f in fs if os.path.exists(f)]
            if fs:
                ims = [Image.open(f).convert("P", palette=Image.ADAPTIVE) for f in fs]
                ims[0].save(os.path.join(ANIM, "grain_" + k + ".gif"), save_all=True, append_images=ims[1:], duration=80, loop=0)
    dn = {r["file"] for r in new}
    rows = sorted([r for r in old if r["file"] not in dn] + new, key=lambda r: r["t_ns"])
    def write_tab(fn, tab, cols):
        ks = ("id",) if "id" in cols else ("grain_i", "grain_j")
        tab = sorted(tab, key=lambda r: (float(r["t_ns"]),) + tuple(int(float(r[k])) for k in ks))
        with open(fn, "w", newline="") as f:
            wr = csv.writer(f); wr.writerow(cols)
            for r in tab: wr.writerow([("%.6g" % r[c]) if isinstance(r[c], float) else r[c] for c in cols])
    if gtab:
        write_tab(gtf, gtab, ["t_ns", "step", "file", "id", "plateau_area_um2", "territory_area_um2", "cx_um", "cy_um", "theta_deg",
                              "theta_std_deg", "theta_sym_deg", "orig_label", "n_neighbours"])
    if btab:
        write_tab(btf, btab, ["t_ns", "step", "file", "grain_i", "grain_j", "theta_i_deg", "theta_j_deg", "misorient_deg", "misorient_sym_deg",
                              "length_um", "width_um", "mean_phi_at_edge", "mid_x_um", "mid_y_um"])
    if not rows:
        sys.exit("grain_post: no frame could be analysed")
    cols = ["t_ns", "step", "file"] + [k for k in rows[0] if k not in ("t_ns", "step", "file")]
    with open(hf, "w", newline="") as f:
        wr = csv.writer(f); wr.writerow(cols)
        for r in rows: wr.writerow([("%.6g" % r[c]) if isinstance(r[c], float) else r[c] for c in cols])
    # detail of the first and last frame (redone only if that frame is new or its files are missing)
    st = json.load(open(sf)) if os.path.exists(sf) else {}
    det = st.get("detail", {}) if st.get("params") == par else {}
    for name, (t, f) in (("first", FR[0]), ("last", FR[-1])):
        pfx = os.path.join(OUT, "grains_" + name)
        if det.get(name) != os.path.basename(f) or not os.path.exists(pfx + "_misorientation.png"):
            try:
                detail(f, pfx); det[name] = os.path.basename(f)
            except Exception as e:
                print("grain_post: detail of the %s frame failed: %s" % (name, e))
    json.dump(dict(params=par, frames=[r["file"] for r in rows], detail=det, lim=lim), open(sf, "w"))
    plot(rows, os.path.join(OUT, "grain_history.png"))
    r = rows[-1]; names = [n for n, _, _ in gc.bins_of(PAR)]
    print("grain_post: last frame t = %.3f ns: %d grains (%d >= 1 um^2), boundary length %.1f um; raw %%: %s   [%.0f s]"
          % (r["t_ns"], r["n_grains"], r["n_ge1um2"], r["boundary_length_um"],
             "  ".join("%s %.1f" % (n, r["pct_raw_" + n]) for n in names), time.time() - t0))
