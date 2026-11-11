#!/usr/bin/env python3
"""Grain-evolution summary figures from grain_post.py outputs (deliver_<tag>/grains/).

Per run (grain_history.csv + grain_table_all.csv):
  refinement factor  R(t)  = d0 / d(t),  d = mean equivalent-circle diameter of grain territories
                     (also area-based  A0/A(t) = N(t)/N0 is the count increase)
  grain count        N(t) and N(t)/N0
  mean grain size    mean +- std of equivalent diameter d = sqrt(4A/pi) over grains at each time
                     (--area territory [default, tiles the domain] | plateau [resolved interior only])
  misorientation     boundary-LENGTH fractions in VL [2,5), L [5,15), H >=15 deg  (and <2 deg, grey)
                     (--misor sym [default, cubic/90-deg reduced] | raw)

Usage
  python3 grain_summary.py deliver_A/grains deliver_B/grains ... --out summary   # any list of grain folders
  python3 grain_summary.py --root .. --angles 0 45 90 --ups 25 50 100 150 250 400 --out ../summary_ladder
Options: --prefix deliver --suffix _<RUNTAG>   (folders <prefix>_a<A>_up<U>ms<suffix>/grains)   --smooth N (running mean, frames)
         --t0 (reference time for N0,d0; default first frame)
Main output: end_grid.png = angle x velocity grids of end-of-run values (--tavg N averages the last N ns).
Outputs: per-run  <out>/<tag>_summary.png ; with --lines  <out>/compare_a<ANGLE>_{metrics,misorient}.png
         <out>/summary_table.csv  (final / peak values per run)
"""
import argparse, glob, os, re, sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PFX = "deliver"
BINS = [("<2", "<2°", "#9a9a9a"), ("[2,5)", "VL 2–5°", "#4c78a8"),
        ("[5,15)", "L 5–15°", "#f2a541"), (">=15", "H >15°", "#c0392b")]


def smooth(y, n):
    y = np.asarray(y, float)
    if n <= 1:
        return y
    k = np.ones(n) / n
    p = np.pad(y, (n // 2, n - 1 - n // 2), mode="edge")
    return np.convolve(p, k, mode="valid")


def load(d, area="territory", misor="sym"):
    h = pd.read_csv(os.path.join(d, "grain_history.csv")).sort_values("t_ns").reset_index(drop=True)
    g = pd.read_csv(os.path.join(d, "grain_table_all.csv"))
    col = "territory_area_um2" if area == "territory" else "plateau_area_um2"
    rows = []
    for t, s in g.groupby("t_ns"):
        a = s[col].to_numpy(float)
        a = a[a > 0]
        dd = np.sqrt(4 * a / np.pi) if len(a) else np.array([np.nan])
        rows.append((t, len(a), dd.mean(), dd.std(ddof=0), a.mean() if len(a) else np.nan))
    sz = pd.DataFrame(rows, columns=["t_ns", "n", "d_mean", "d_std", "a_mean"])
    h = h.merge(sz, on="t_ns", how="left")
    p = "pct_sym_" if misor == "sym" else "pct_raw_"
    for k, _, _ in BINS:
        h["m" + k] = h[p + k]
    return h


def derived(h, t0=None):
    i0 = 0 if t0 is None else int((h.t_ns - t0).abs().idxmin())
    N0 = h.n_grains.iloc[i0]
    d0 = h.d_mean.iloc[i0]
    h = h.copy()
    h["N_ratio"] = h.n_grains / N0
    h["refine"] = d0 / h.d_mean
    h.attrs.update(N0=N0, d0=d0)
    return h


def run_fig(h, tag, out, sm, area, misor):
    t = h.t_ns.to_numpy()
    fig, ax = plt.subplots(2, 2, figsize=(11, 7.2), sharex=True)
    a = ax[0, 0]
    a.plot(t, smooth(h.refine, sm), color="#1b6ca8", lw=1.6)
    a.axhline(1, color="k", lw=.6, ls=":")
    a.set(ylabel=r"refinement factor  $d_0/d(t)$", title="Grain refinement factor")
    a = ax[0, 1]
    a.plot(t, h.n_grains, color="#2e8b57", lw=1.6, label="N(t)")
    a.set(ylabel="grains N", title="Grain count (N/N$_0$ on right)")
    b = a.twinx()
    b.plot(t, h.N_ratio, color="none")
    b.set_ylim(np.array(a.get_ylim()) / h.attrs["N0"])
    b.set_ylabel(r"$N/N_0$")
    a = ax[1, 0]
    m, s = smooth(h.d_mean, sm), smooth(h.d_std, sm)
    a.plot(t, m, color="#6a3d9a", lw=1.6, label="mean")
    a.fill_between(t, np.maximum(m - s, 0), m + s, color="#6a3d9a", alpha=.2, label="±1 std")
    a.set(ylabel=r"equivalent diameter $d$ [µm]", xlabel="t [ns]",
          title=f"Mean grain size ± std ({area} area)")
    a.legend(frameon=False)
    a = ax[1, 1]
    ys = [smooth(h["m" + k], sm) for k, _, _ in BINS]
    a.stackplot(t, ys, colors=[c for _, _, c in BINS], labels=[l for _, l, _ in BINS], alpha=.9)
    a.set(ylim=(0, 100), ylabel="boundary length [%]", xlabel="t [ns]",
          title=f"Misorientation distribution ({misor})")
    a.legend(frameon=False, loc="upper right", fontsize=8, ncol=2)
    for x in ax.ravel():
        x.grid(alpha=.25)
    fig.suptitle(tag, fontsize=10)
    fig.tight_layout()
    fig.savefig(out, dpi=170)
    plt.close(fig)


def parse(tag):
    m = re.search(r"(?:^|_)a(\d+)_up(\d+)ms", tag)
    return (int(m.group(1)), int(m.group(2))) if m else (None, None)


def compare_figs(runs, ang, outdir, sm, area, misor):
    runs = sorted(runs, key=lambda r: r[2])
    cm = plt.get_cmap("viridis")
    cols = {r[2]: cm(i / max(len(runs) - 1, 1) * .9) for i, r in enumerate(runs)}
    f, ax = plt.subplots(1, 3, figsize=(14, 4.2), sharex=True)
    spec = [("refine", r"refinement factor $d_0/d$"), ("n_grains", "grain count N"),
            ("d_mean", r"mean diameter $d$ [µm] (±std band)")]
    for (up, h) in [(r[2], r[1]) for r in runs]:
        t = h.t_ns
        ax[0].plot(t, smooth(h.refine, sm), color=cols[up], label=f"{up} m/s")
        ax[1].plot(t, h.n_grains, color=cols[up])
        m, s = smooth(h.d_mean, sm), smooth(h.d_std, sm)
        ax[2].plot(t, m, color=cols[up])
        ax[2].fill_between(t, np.maximum(m - s, 0), m + s, color=cols[up], alpha=.12)
    for x, (_, yl) in zip(ax, spec):
        x.set(xlabel="t [ns]", ylabel=yl)
        x.grid(alpha=.25)
    ax[0].axhline(1, color="k", lw=.6, ls=":")
    ax[0].legend(frameon=False, title=f"angle {ang}°")
    f.tight_layout()
    f.savefig(os.path.join(outdir, f"compare_a{ang}_metrics.png"), dpi=170)
    plt.close(f)
    f, ax = plt.subplots(1, 3, figsize=(14, 4.2), sharex=True, sharey=True)
    for x, (k, l, _) in zip(ax, BINS[1:]):
        for r in runs:
            x.plot(r[1].t_ns, smooth(r[1]["m" + k], sm), color=cols[r[2]], label=f"{r[2]} m/s")
        x.set(title=f"{l}  ({misor})", xlabel="t [ns]")
        x.grid(alpha=.25)
    ax[0].set_ylabel("boundary length [%]")
    ax[0].legend(frameon=False, title=f"angle {ang}°")
    f.tight_layout()
    f.savefig(os.path.join(outdir, f"compare_a{ang}_misorient.png"), dpi=170)
    plt.close(f)


def end_value(h, tavg):
    """end-of-run values: last frame, or mean over the last tavg ns"""
    w = h[h.t_ns >= h.t_ns.max() - tavg] if tavg > 0 else h.iloc[[-1]]
    return w.mean(numeric_only=True)


def grid_fig(runs, out, tavg, ups, angs, misor, area, size=(12.8, 7.2), suptitle=False):
    """angle x velocity grids of end-of-run values, numbers printed in each cell"""
    E = {(r[0], r[2]): end_value(r[1], tavg) for r in runs}
    N0 = next(iter(runs))[1].attrs["N0"]
    panels = [("refine", r"refinement factor $d_0/d$", "Blues", "{:.2f}"),
              ("n_grains", f"grain count N  (N$_0$={N0})", "Greens", "{:.0f}"),
              ("d_mean", r"mean diameter $d$ ± std [µm]", "Purples_r", None),
              ("m[2,5)", f"VL 2–5° [% length] ({misor})", "Blues", "{:.0f}"),
              ("m[5,15)", "L 5–15° [% length]", "Oranges", "{:.0f}"),
              ("m>=15", "H >15° [% length]", "Reds", "{:.0f}")]
    nc = 3
    nr = -(-len(panels) // nc)
    f, axs = plt.subplots(nr, nc, figsize=tuple(size))
    axs = np.atleast_1d(axs).ravel()
    init = {}
    for r in runs:
        init.setdefault(r[0], r[1].iloc[0])
    for ax, (k, ptitle, cmap, fmt) in zip(axs, panels):
        off = 0
        M = np.full((len(angs), len(ups) + off), np.nan)
        S = np.full_like(M, np.nan)
        for i, A in enumerate(angs):
            if off and init:                       # t=0 state is the same microstructure for all runs; use own angle if present
                r0 = init.get(A, next(iter(init.values())))
                M[i, 0], S[i, 0] = r0["d_mean"], r0["d_std"]
            for j, U in enumerate(ups):
                if (A, U) in E:
                    M[i, j + off] = E[(A, U)][k]
                    if k == "d_mean":
                        S[i, j + off] = E[(A, U)]["d_std"]
        im = ax.imshow(np.ma.masked_invalid(M), cmap=cmap, aspect="auto")
        lo, hi = np.nanmin(M), np.nanmax(M)
        im.set_clim(lo, hi if hi > lo else lo + 1)
        for i in range(M.shape[0]):
            for j in range(M.shape[1]):
                if np.isnan(M[i, j]):
                    ax.text(j, i, "–", ha="center", va="center", color="#999")
                    continue
                txt = f"{M[i, j]:.2f}\n±{S[i, j]:.2f}" if k == "d_mean" else fmt.format(M[i, j])
                r, g_, b_, _ = im.cmap(im.norm(M[i, j]))
                dark = (0.299 * r + 0.587 * g_ + 0.114 * b_) < .5
                ax.text(j, i, txt, ha="center", va="center", fontsize=12,
                        color="white" if dark else "black")
        if off:
            ax.axvline(.5, color="k", lw=2)
        ax.set_xticks(range(len(ups) + off), (["t=0"] if off else []) + [str(u) for u in ups])
        ax.set_yticks(range(len(angs)), [f"{a}°" for a in angs])
        ax.set_xlabel("plate speed [m/s]" + ("  (end of run)" if off else ""), fontsize=11)
        ax.set_ylabel("angle", fontsize=11)
        ax.tick_params(labelsize=11)
        if k == "d_mean" and init:
            r0 = next(iter(init.values()))          # t=0 microstructure is common to all runs
            ptitle += f"\ninitial (t=0): {r0['d_mean']:.2f} ± {r0['d_std']:.2f} µm"
        ax.set_title(ptitle, fontsize=13)
        ax.tick_params(length=0)
    for ax in axs[len(panels):]:
        ax.axis("off")
    when = f"mean of last {tavg:g} ns" if tavg > 0 else "last frame"
    if suptitle:
        f.suptitle(f"Grain evolution, end of run ({when}; size from {area} area)", fontsize=11)
    f.tight_layout()
    f.savefig(out, dpi=170)
    f.savefig(os.path.splitext(out)[0] + ".pdf")
    plt.close(f)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("dirs", nargs="*")
    ap.add_argument("--root", default=None)
    ap.add_argument("--prefix", default=PFX)
    ap.add_argument("--suffix", default="", help="tag suffix after _up<U>ms (e.g. _gm for RUNTAG=gm, _serial)")
    ap.add_argument("--angles", type=int, nargs="*", default=[0, 45, 90])
    ap.add_argument("--ups", type=int, nargs="*", default=[25, 50, 100, 150, 250, 400])
    ap.add_argument("--out", default="grain_summary")
    ap.add_argument("--smooth", type=int, default=1)
    ap.add_argument("--area", choices=["territory", "plateau"], default="territory")
    ap.add_argument("--misor", choices=["sym", "raw"], default="sym")
    ap.add_argument("--tavg", type=float, default=0, help="average end values over the last TAVG ns (0 = last frame)")
    ap.add_argument("--lines", action="store_true", help="also write the per-angle time-series overlays")
    ap.add_argument("--no-per-run", action="store_true")
    ap.add_argument("--size", type=float, nargs=2, default=[12.8, 7.2], help="end_grid figure size in inches (default 16:9 slide; 4:3 slide: 10 7.5)")
    ap.add_argument("--title", action="store_true", help="add a suptitle to end_grid")
    ap.add_argument("--t0", type=float, default=None)
    a = ap.parse_args()
    dirs = list(a.dirs)
    if a.root:
        for A in a.angles:
            for U in a.ups:
                base = os.path.join(a.root, f"{a.prefix}_a{A}_up{U}ms{a.suffix}")
                for d in (os.path.join(base, "grains"), base):   # deliver_<tag>/grains, or an older grains_<tag>/ folder
                    if os.path.exists(os.path.join(d, "grain_history.csv")):
                        dirs.append(d); break
                else:
                    print(f"[skip] missing {base}/grains", file=sys.stderr)
    if not dirs:
        sys.exit("no run folders found")
    os.makedirs(a.out, exist_ok=True)
    runs, tab = [], []
    for d in dirs:
        nd = os.path.normpath(d)
        tag = os.path.basename(nd)
        if tag == "grains": tag = os.path.basename(os.path.dirname(os.path.abspath(nd)))   # deliver_<tag>/grains
        tag = re.sub(r"^(grains|deliver)_", "", tag)
        try:
            h = derived(load(d, a.area, a.misor), a.t0)
        except Exception as e:
            print(f"[skip] {d}: {e}", file=sys.stderr)
            continue
        ang, up = parse(tag)
        if not a.no_per_run:
            run_fig(h, tag, os.path.join(a.out, f"{tag}_summary.png"), a.smooth, a.area, a.misor)
        runs.append((ang, h, up))
        L = h.iloc[-1]
        tab.append(dict(tag=tag, angle=ang, up_ms=up, t_end=L.t_ns, N0=h.attrs["N0"], N_end=L.n_grains,
                        N_peak=h.n_grains.max(), N_ratio_end=L.N_ratio, d0_um=h.attrs["d0"],
                        d_end_um=L.d_mean, d_std_end_um=L.d_std, refine_end=L.refine,
                        refine_max=h.refine.max(), t_refine_max=h.t_ns[h.refine.idxmax()],
                        VL_end=L["m[2,5)"], L_end=L["m[5,15)"], H_end=L["m>=15"], lt2_end=L["m<2"]))
    if not runs:
        sys.exit("nothing plotted")
    pd.DataFrame(tab).to_csv(os.path.join(a.out, "summary_table.csv"), index=False, float_format="%.4g")
    byang = {}
    for ang, h, up in runs:
        if ang is not None:
            byang.setdefault(ang, []).append((ang, h, up))
    angs, ups = list(a.angles), list(a.ups)   # full grid; runs not (yet) available show as blank
    grid_fig([r for r in runs if r[0] is not None], os.path.join(a.out, "end_grid.png"), a.tavg, ups, angs, a.misor, a.area, a.size, a.title)
    for ang, rr in (byang.items() if a.lines else []):
        if len(rr) > 1:
            compare_figs(rr, ang, a.out, a.smooth, a.area, a.misor)
    print(f"wrote {len(runs)} run figure(s) to {a.out}")


if __name__ == "__main__":
    main()
