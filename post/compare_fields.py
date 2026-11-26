#!/usr/bin/env python3
"""Compare the plastic fields of two runs frame by frame, e.g. inertial (j = 0.01) against inertia-less (j = 0).

  python3 post/compare_fields.py RUN_A RUN_B [--out DIR] [--fields f1 f2 ...] [--every N] [--tmax T]
                                 [--maps t1 t2 ...] [--mapfields f1 f2 ...] [--labels A B] [--increments]

RUN_A is the reference (e.g. out_a45_up400ms_v2_kt3b), RUN_B the other run (e.g. out_a45_up400ms_j_kt3b). Both must use the same
mesh (the VTU cells are compared one by one, in the reference configuration). Frames are matched by time from the .pvd files.
The two set-ups relax separately, so their t = 0 frames can already differ; --increments compares the change since t = 0 instead.

Output in --out (default cmp_<A>_vs_<B>):
  compare.csv          one row per (time, field): area-weighted mean in A and B and their difference in % of |mean A|; RMS of B - A
                       absolute and in % of RMS(A); max |B - A| and its 99th percentile, absolute and in % of max |A|; the
                       fraction of the area where |B - A| > 10 % of max |A|
  compare_means.png    area-weighted mean of every field over time, A and B
  compare_means_pct.png  difference of the means, (mean B - mean A)/|mean A| in %, over time
  compare_diff.png     RMS(B - A)/RMS(A) in % over time, one curve per field (log scale)
  maps_t<t>.png        A, B, B - A in % of max |A|, and the local difference (B - A)/|A| in % where |A| > 5 % of max |A|,
                       at the requested times (default: four times spread over the common range)
Needs numpy, meshio, matplotlib.
"""
import argparse, csv, glob, os, re, sys
import numpy as np

try:
    import meshio
except ImportError:
    sys.exit("compare_fields.py needs meshio (pip install --user meshio)")

PLASTIC = ["slip_path_phys", "slip_path_eff", "p1", "p2", "p3", "wstar", "rhoS", "rhoG",
           "Fp_11", "Fp_12", "Fp_21", "Fp_22", "Fp_dev", "D_rate"]
MAPS = ["slip_path_phys", "rhoS", "wstar", "Fp_dev"]


def read_pvd(d):
    out = {}
    for pvd in glob.glob(os.path.join(d, "*.pvd")):
        for tm, fn in re.findall(r'timestep="([^"]+)"\s+file="([^"]+)"', open(pvd).read()):
            out[round(float(tm), 6)] = os.path.join(d, os.path.basename(fn))
    if not out:
        sys.exit(f"no .pvd with frames in {d}")
    return out


def load(path, fields):
    m = meshio.read(path)
    k = [i for i, c in enumerate(m.cells) if c.type == "triangle"]
    if not k:
        sys.exit(f"no triangles in {path}")
    k = k[0]
    tri = m.cells[k].data
    f = {}
    for name in fields:
        if name == "Fp_dev":   # |F_p - I| (Frobenius): total plastic distortion
            a = [np.asarray(m.cell_data[c][k], float).ravel() for c in ("Fp_11", "Fp_12", "Fp_21", "Fp_22")]
            f[name] = np.sqrt((a[0] - 1) ** 2 + a[1] ** 2 + a[2] ** 2 + (a[3] - 1) ** 2)
        elif name in m.cell_data:
            f[name] = np.asarray(m.cell_data[name][k], float).ravel()
        elif name in m.point_data:   # point field -> element mean
            f[name] = np.asarray(m.point_data[name], float).reshape(len(m.points), -1)[:, 0][tri].mean(axis=1)
    return m.points[:, :2], tri, f


def areas(p, t):
    a, b, c = p[t[:, 0]], p[t[:, 1]], p[t[:, 2]]
    return 0.5 * np.abs((b[:, 0] - a[:, 0]) * (c[:, 1] - a[:, 1]) - (b[:, 1] - a[:, 1]) * (c[:, 0] - a[:, 0]))


def wq(x, w, q):
    o = np.argsort(x); c = np.cumsum(w[o]) - 0.5 * w[o]; c /= w.sum()
    return float(np.interp(q, c, x[o]))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_a"); ap.add_argument("run_b")
    ap.add_argument("--out"); ap.add_argument("--fields", nargs="+", default=PLASTIC)
    ap.add_argument("--every", type=int, default=1, help="use every N-th common frame")
    ap.add_argument("--tmax", type=float, default=np.inf)
    ap.add_argument("--maps", nargs="*", type=float, help="times (ns) for field maps; default: four spread over the run")
    ap.add_argument("--mapfields", nargs="+", default=MAPS)
    ap.add_argument("--labels", nargs=2, help="legend labels for A and B (default: the folder names)")
    ap.add_argument("--increments", action="store_true",
                    help="compare the change of every field since the first common frame (removes a difference in the relaxed states)")
    a = ap.parse_args()
    na, nb = (os.path.basename(os.path.normpath(r)) for r in (a.run_a, a.run_b))
    la, lb = a.labels or (na, nb)
    out = a.out or f"cmp_{na}_vs_{nb}"
    os.makedirs(out, exist_ok=True)

    fa, fb = read_pvd(a.run_a), read_pvd(a.run_b)
    times = sorted(t for t in set(fa) & set(fb) if t <= a.tmax)[::a.every]
    if not times:
        sys.exit("no common frame times")
    print(f"{len(times)} common frames, t = {times[0]} .. {times[-1]} ns")
    fields = list(dict.fromkeys(a.fields + [f for f in a.mapfields if f not in a.fields]))

    rows, series, w = [], {f: {k: [] for k in ("t", "ma", "mb", "rel", "mp")} for f in a.fields}, None
    maps = a.maps if a.maps is not None else [times[int(round(i * (len(times) - 1) / 3))] for i in range(4)]
    maps = sorted({min(times, key=lambda t: abs(t - m)) for m in maps})
    snaps = {}
    for t in times:
        p, tri, A = load(fa[t], fields)
        p2, tri2, B = load(fb[t], fields)
        if len(tri) != len(tri2) or not np.array_equal(tri, tri2):
            sys.exit(f"the meshes of {fa[t]} and {fb[t]} differ: fields cannot be compared cell by cell")
        if w is None:
            w = areas(p, tri)
            A0, B0 = ({f: v.copy() for f, v in A.items()}, {f: v.copy() for f, v in B.items()})
        if a.increments:
            A = {f: v - A0[f] for f, v in A.items()}; B = {f: v - B0[f] for f, v in B.items()}
        W = w.sum()
        for f in a.fields:
            if f not in A or f not in B:
                continue
            x, y = A[f], B[f]; d = y - x
            rmsA = np.sqrt(np.sum(w * x * x) / W); rmsd = np.sqrt(np.sum(w * d * d) / W)
            mxA = np.abs(x).max()
            mA, mB = np.sum(w * x) / W, np.sum(w * y) / W
            mx, p99 = np.abs(d).max(), wq(np.abs(d), w, 0.99)
            pct = lambda v, ref: 100.0 * v / ref if ref > 0 else np.nan
            r = dict(t_ns=t, field=f, mean_A=mA, mean_B=mB, mean_diff_pct=pct(mB - mA, abs(mA)),
                     rms_diff=rmsd, rms_diff_pct=pct(rmsd, rmsA),
                     max_abs_diff=mx, max_abs_diff_pct_of_max=pct(mx, mxA), p99_abs_diff=p99, p99_abs_diff_pct_of_max=pct(p99, mxA),
                     area_frac_10pct=(np.sum(w[np.abs(d) > 0.1 * mxA]) / W if mxA > 0 else np.nan))
            rows.append(r)
            s = series[f]; s["t"].append(t); s["ma"].append(r["mean_A"]); s["mb"].append(r["mean_B"]); s["rel"].append(r["rms_diff_pct"]); s["mp"].append(r["mean_diff_pct"])
        if t in maps:
            snaps[t] = (p, tri, A, B)
        print(f"  t = {t:8.3f} ns done", end="\r")
    print()

    with open(os.path.join(out, "compare.csv"), "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(rows[0])); wr.writeheader(); wr.writerows(rows)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.tri as mtri
    fl = [f for f in a.fields if series[f]["t"]]
    nc = 4; nr = int(np.ceil(len(fl) / nc))
    fig, ax = plt.subplots(nr, nc, figsize=(4 * nc, 2.8 * nr), squeeze=False)
    for i, f in enumerate(fl):
        s = series[f]; x = ax[i // nc][i % nc]
        x.plot(s["t"], s["ma"], label=la, lw=1.5); x.plot(s["t"], s["mb"], "--", label=lb, lw=1.5)
        x.set_title(f, fontsize=9); x.set_xlabel("t (ns)", fontsize=8); x.tick_params(labelsize=7); x.grid(alpha=0.3)
    for j in range(len(fl), nr * nc):
        ax[j // nc][j % nc].axis("off")
    ax[0][0].legend(fontsize=7)
    fig.suptitle("area-weighted mean of the plastic fields", fontsize=10); fig.tight_layout()
    fig.savefig(os.path.join(out, "compare_means.png"), dpi=150); plt.close(fig)

    fig, x = plt.subplots(figsize=(8, 4.8))
    for f in fl:
        x.plot(series[f]["t"], series[f]["mp"], label=f, lw=1.3)
    x.axhline(0, color="k", lw=0.6); x.set_xlabel("t (ns)"); x.set_ylabel("(mean B − mean A) / |mean A|  (%)"); x.grid(alpha=0.3)
    x.set_title(f"difference of the area-weighted means, A = {la}, B = {lb}", fontsize=10); x.legend(fontsize=7, ncol=2)
    fig.tight_layout(); fig.savefig(os.path.join(out, "compare_means_pct.png"), dpi=150); plt.close(fig)

    fig, x = plt.subplots(figsize=(8, 4.8))
    for f in fl:
        s = series[f]; r = np.array(s["rel"], float)
        if np.any(np.isfinite(r) & (r > 0)):
            x.semilogy(s["t"], r, label=f, lw=1.3)
    x.set_xlabel("t (ns)"); x.set_ylabel("RMS(B − A) / RMS(A)  (%)"); x.grid(alpha=0.3, which="both")
    x.set_title(f"relative difference, A = {la}, B = {lb}", fontsize=10); x.legend(fontsize=7, ncol=2)
    fig.tight_layout(); fig.savefig(os.path.join(out, "compare_diff.png"), dpi=150); plt.close(fig)

    for t, (p, tri, A, B) in snaps.items():
        mf = [f for f in a.mapfields if f in A and f in B]
        tr = mtri.Triangulation(p[:, 0], p[:, 1], tri)
        fig, ax = plt.subplots(len(mf), 4, figsize=(16, 3.6 * len(mf)), squeeze=False)
        for i, f in enumerate(mf):
            lo, hi = min(A[f].min(), B[f].min()), max(A[f].max(), B[f].max())
            for j, (v, ttl) in enumerate(((A[f], la), (B[f], lb))):
                im = ax[i][j].tripcolor(tr, facecolors=v, vmin=lo, vmax=hi, cmap="viridis")
                ax[i][j].set_title(f"{f}: {ttl}", fontsize=9); fig.colorbar(im, ax=ax[i][j], shrink=0.8)
            d = B[f] - A[f]; mA = np.abs(A[f]).max() or 1.0
            dp = 100.0 * d / mA; m = np.abs(dp).max() or 1.0
            im = ax[i][2].tripcolor(tr, facecolors=dp, vmin=-m, vmax=m, cmap="RdBu_r")
            ax[i][2].set_title(f"{f}: B − A, % of max |A|", fontsize=9); fig.colorbar(im, ax=ax[i][2], shrink=0.8)
            loc = np.ma.masked_where(np.abs(A[f]) <= 0.05 * mA, 100.0 * d / np.where(np.abs(A[f]) > 0, np.abs(A[f]), 1.0))
            cmap = matplotlib.colormaps["RdBu_r"].copy(); cmap.set_bad("#d9d9d9")
            im = ax[i][3].tripcolor(tr, facecolors=np.ma.clip(loc, -100, 100), vmin=-100, vmax=100, cmap=cmap)
            ax[i][3].set_title(f"{f}: (B − A)/|A|, %  (grey: |A| ≤ 5 % of max)", fontsize=9); fig.colorbar(im, ax=ax[i][3], shrink=0.8)
            for j in range(4):
                ax[i][j].set_aspect("equal"); ax[i][j].set_xticks([]); ax[i][j].set_yticks([])
        fig.suptitle(f"t = {t:g} ns   (A = {la}, B = {lb})", fontsize=10); fig.tight_layout()
        fig.savefig(os.path.join(out, f"maps_t{t:g}.png"), dpi=130); plt.close(fig)
    print(f"wrote {out}/compare.csv, compare_means.png, compare_means_pct.png, compare_diff.png and {len(snaps)} map figures")


if __name__ == "__main__":
    main()
