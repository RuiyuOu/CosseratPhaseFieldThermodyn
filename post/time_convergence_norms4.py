#!/usr/bin/env python3
"""time_convergence_norms4.py -- time-step convergence with FOUR runs: dt, dt/2, dt/4 compared against dt/8 as the reference.

usage:  python3 time_convergence_norms4.py --runs out_<dt> out_<dt/2> out_<dt/4> out_<dt/8> [--out post/time_convergence_dt8] [--every 0.5] [--np 4]
(any order; sorted by dt, the smallest dt is the reference).  Companion of time_convergence_norms.py: same frame reader, area weights,
norms (L2, P99, Linf), the same two denominators ("change" = ||f_ref - f_0||, "magnitude" = ||f_ref||) and the same space-time norms.

For every field and norm:  e_k = ||f_k - f_ref|| / denominator for k = dt, dt/2, dt/4  (three data points).
Observed order from the data points, two versions:
  p_naive(a,b) = log2(e_a / e_b)                  slope of the error against the reference, biased because the reference itself has an error
  p_corr       solves e_a/e_b = (h_a^p - h_r^p)/(h_b^p - h_r^p) exactly for the step ratios used (h_r = dt/8):
                 p_corr(dt/2,dt/4) = log2(e_2/e_3 - 1),   p_corr(dt,dt/2) by bisection on (8^p-1)/(4^p-1) = e_1/e_2.
e_true_dt = e_1 + e_3/(2^p - 1)  estimates the error of the production dt including the error of the reference itself (p = p_corr(dt/2,dt/4)
clipped to [0.5, 2]).  Without a relaxed frame in the run folders the reference state f_0 is the earliest frame common to all four runs.

Outputs in --out:  norm_errors4.csv (per time), norm_summary4.csv (space-time), vtu_means4.csv, norm_summary4.json, tc_norms4.png
"""
import argparse, json, os, sys
from multiprocessing import Pool
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import time_convergence as tc
import time_convergence_norms as tn

NORMS, MEANS = tn.NORMS, tn.MEANS


def p_corr_23(e2, e3):
    return float(np.log2(e2 / e3 - 1.0)) if np.isfinite(e2) and np.isfinite(e3) and e3 > 0 and e2 / e3 > 1.0 else float("nan")


def p_corr_12(e1, e2):
    """solve (8^p-1)/(4^p-1) = e1/e2 for p in (0.05, 6]; the left side rises monotonically from 1.5 (p -> 0) to ~2^p"""
    if not (np.isfinite(e1) and np.isfinite(e2)) or e2 <= 0: return float("nan")
    r = e1 / e2; f = lambda p: (8.0 ** p - 1.0) / (4.0 ** p - 1.0) - r
    lo, hi = 0.05, 6.0
    if f(lo) * f(hi) > 0: return float("nan")
    for _ in range(80):
        m = 0.5 * (lo + hi)
        if f(lo) * f(m) <= 0: hi = m
        else: lo = m
    return float(0.5 * (lo + hi))


def p_naive(ea, eb):
    return float(np.log2(ea / eb)) if np.isfinite(ea) and np.isfinite(eb) and ea > 0 and eb > 0 else float("nan")


def job(args):
    t, fs, f0 = args                                          # fs = frames of dt, dt/2, dt/4, dt/8 (reference)
    P, T, F1 = tc.load_fields(fs[0]); _, _, F2 = tc.load_fields(fs[1]); _, _, F3 = tc.load_fields(fs[2]); _, _, FR = tc.load_fields(fs[3])
    _, _, F0 = tc.load_fields(f0)
    wn, wc = tc.weights(P, T)
    rows, means = [], []
    for k in tc.NODAL + tc.CELL:
        w = wn if k in tc.NODAL else wc
        rc, rm = tn.norms(FR[k] - F0[k], w), tn.norms(FR[k], w)
        d = [tn.norms(F[k] - FR[k], w) for F in (F1, F2, F3)]
        for n in NORMS:
            rows.append(dict(t_ns=t, field=k, norm=n, ref_change=rc[n], ref_mag=rm[n], d1=d[0][n], d2=d[1][n], d3=d[2][n]))
        if k in MEANS:
            means.append(dict(t_ns=t, field=k, **{"m%d" % (i + 1): float((w * F[k]).sum() / w.sum()) for i, F in enumerate((F1, F2, F3, FR))}))
    return rows, means


def summarize(d1, d2, d3, rc, rm):
    sc = rc if rc > 0 else float("nan"); sm = rm if rm > 0 else float("nan")
    e1, e2, e3 = d1 / sc, d2 / sc, d3 / sc
    pc = p_corr_23(e2, e3); pcl = min(2.0, max(0.5, pc)) if np.isfinite(pc) else float("nan")
    et = e1 + e3 / (2.0 ** pcl - 1.0) if np.isfinite(pcl) else float("nan")
    return dict(e_dt=e1, e_dt2=e2, e_dt4=e3, p_naive_12=p_naive(e1, e2), p_naive_23=p_naive(e2, e3),
                p_corr_12=p_corr_12(e1, e2), p_corr_23=pc, e_true_dt=et, em_dt=d1 / sm, em_dt2=d2 / sm, em_dt4=d3 / sm)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--runs", nargs=4, required=True); ap.add_argument("--out", default="post/time_convergence_dt8")
    ap.add_argument("--every", type=float, default=0.5); ap.add_argument("--np", type=int, default=int(os.environ.get("NP", 4)))
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
    runs = sorted(a.runs, key=tc.dt_of, reverse=True); dts = [tc.dt_of(d) for d in runs]
    print("runs:", ", ".join("%s (dt=%g)" % (tc.tag_of(d), x) for d, x in zip(runs, dts)), " -> reference = last")
    ratios = [dts[i] / dts[i + 1] for i in range(3)]
    if any(abs(r - 2.0) > 1e-6 for r in ratios): sys.exit("the four steps must halve each time (ratios %s)" % ratios)
    Fr = [tc.frames(d) for d in runs]
    tf = sorted(set.intersection(*[set(f) for f in Fr]))
    have0 = all(0.0 in f for f in Fr)
    if have0: t0 = 0.0
    else:
        c = [t for t in tf if t > 0]
        if len(c) < 2: sys.exit("no relaxed frame and fewer than 2 common frames; frames per run: %s" % [len(f) for f in Fr])
        t0 = c[0]
        print("note: no relaxed frame in the run folders; reference state = the reference-run frame at t = %g ns (denominators: change since then; errors start after it)" % t0)
    pick = [t for t in tf if t > t0 and abs(t / a.every - round(t / a.every)) < 1e-6]
    if tf and tf[-1] > t0 and tf[-1] not in pick: pick.append(tf[-1])
    f0 = Fr[3][t0]
    jobs = [(t, [f[t] for f in Fr], f0) for t in pick]
    if have0: jobs.insert(0, (0.0, [f[0.0] for f in Fr], f0))
    with Pool(a.np) as pool: res = pool.map(job, jobs)
    raw = [r for rr, _ in res for r in rr]; means = [r for _, mm in res for r in mm]
    print("frames compared:", len(jobs), "(t = %s ... %s ns)" % (jobs[0][0], jobs[-1][0]))
    rows = []
    for r in raw:
        e = summarize(r["d1"], r["d2"], r["d3"], r["ref_change"], r["ref_mag"]) if r["t_ns"] > t0 else {}
        rows.append({**r, **e})
    keys = sorted({k for r in rows for k in r}, key=lambda k: list(rows[0]).index(k) if k in rows[0] else 99)
    tc.wcsv(os.path.join(a.out, "norm_errors4.csv"), [{k: r.get(k, "") for k in keys} for r in rows])
    summ = []
    for k in tc.NODAL + tc.CELL:
        for n in NORMS:
            rr = sorted([r for r in raw if r["field"] == k and r["norm"] == n], key=lambda r: r["t_ns"])
            if len(rr) < 2: continue
            t = np.array([r["t_ns"] for r in rr]); g = lambda key: np.array([r[key] for r in rr])
            if n == "L2": agg = lambda x: float(np.sqrt(tc.trap(x ** 2, t))); lab = "L2(t;L2)"
            else: agg = lambda x: float(np.max(x)); lab = "Linf(t;%s)" % n
            summ.append(dict(field=k, norm=lab, **summarize(agg(g("d1")), agg(g("d2")), agg(g("d3")), agg(g("ref_change")), agg(g("ref_mag")))))
    tc.wcsv(os.path.join(a.out, "norm_summary4.csv"), summ)
    mrows = []
    for k in MEANS:
        mm = sorted([m for m in means if m["field"] == k], key=lambda m: m["t_ns"])
        if len(mm) < 2: continue
        t = np.array([m["t_ns"] for m in mm]); A = [np.array([m["m%d" % i] for m in mm]) for i in (1, 2, 3, 4)]
        nrm = lambda x: float(np.sqrt(tc.trap(x ** 2, t))); n4 = nrm(A[3])
        mrows.append(dict(quantity="mean_" + k, norm_ref=n4, **summarize(nrm(A[0] - A[3]), nrm(A[1] - A[3]), nrm(A[2] - A[3]), n4, n4),
                          end_dt=A[0][-1], end_dt2=A[1][-1], end_dt4=A[2][-1], end_dt8=A[3][-1]))
    if mrows: tc.wcsv(os.path.join(a.out, "vtu_means4.csv"), mrows)
    json.dump(dict(runs=[tc.tag_of(d) for d in runs], dt=dts, times=[j[0] for j in jobs], space_time=summ, vtu_means=mrows),
              open(os.path.join(a.out, "norm_summary4.json"), "w"), indent=1, default=float)
    print("\nspace-time errors vs the dt/8 reference, relative to the change since the reference state (dt = %g, dt/2, dt/4):" % dts[0])
    print("  %-18s %-12s %9s %9s %9s  %7s %7s  %7s %7s  %10s" % ("field", "norm", "e(dt)", "e(dt/2)", "e(dt/4)", "p12 n", "p23 n", "p12 c", "p23 c", "e_true(dt)"))
    for s in summ:
        print("  %-18s %-12s %9.2e %9.2e %9.2e  %7.2f %7.2f  %7.2f %7.2f  %10.2e" % (s["field"], s["norm"], s["e_dt"], s["e_dt2"], s["e_dt4"],
              s["p_naive_12"], s["p_naive_23"], s["p_corr_12"], s["p_corr_23"], s["e_true_dt"]))
    print("\ndomain integrals from the VTU (history L2(0,T)):")
    for m in mrows:
        print("  %-18s e %9.2e %9.2e %9.2e  p23c %5.2f  e_true(dt) %9.2e" % (m["quantity"], m["e_dt"], m["e_dt2"], m["e_dt4"], m["p_corr_23"], m["e_true_dt"]))
    # figure: error against dt (log-log), L2(t;L2) relative to the change, with slope-1 and slope-2 guides
    sel = [s for s in summ if s["norm"] == "L2(t;L2)"]
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    h = np.array(dts[:3])
    for i, s in enumerate(sel):
        ax[0].loglog(h, [s["e_dt"], s["e_dt2"], s["e_dt4"]], "o-", color="C%d" % (i % 10), lw=1, label=s["field"])
    if sel:
        e0 = np.nanmedian([s["e_dt4"] for s in sel])
        for pp, ls in ((1, "k:"), (2, "k--")): ax[0].loglog(h, e0 * (h / h[-1]) ** pp, ls, lw=.8, label="slope %d" % pp)
    ax[0].set_xlabel("dt (ns)"); ax[0].set_ylabel("L2(t;L2) error vs dt/8, rel. to change"); ax[0].legend(fontsize=6, ncol=2); ax[0].grid(alpha=.3, which="both")
    for i, s in enumerate(sel):
        ax[1].plot([1, 2], [s["p_corr_12"], s["p_corr_23"]], "o-", color="C%d" % (i % 10), lw=1)
    ax[1].axhspan(1, 2, color="0.9"); ax[1].set_xticks([1, 2]); ax[1].set_xticklabels(["dt -> dt/2", "dt/2 -> dt/4"]); ax[1].set_ylabel("observed order (corrected for the reference error)")
    ax[1].grid(alpha=.3); fig.tight_layout(); fig.savefig(os.path.join(a.out, "tc_norms4.png"), dpi=130); plt.close(fig)
    print("written to", a.out)


if __name__ == "__main__":
    main()
