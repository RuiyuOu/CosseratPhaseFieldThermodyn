#!/usr/bin/env python3
"""time_convergence_global4.py -- global-history time-step convergence with FOUR runs: dt, dt/2, dt/4 against dt/8 as the reference.

usage:  python3 time_convergence_global4.py --runs out_<dt> out_<dt/2> out_<dt/4> out_<dt/8> [--out post/time_convergence_dt8]
(any order; sorted by dt, the smallest dt is the reference).  Companion of time_convergence_norms4.py (same orders and error formulas, imported from it).

For every history in <tag>_log.csv (load_GPa, KE_trans, E_elastic, E_stored, Q_heat, W_ext, D_cum, mean_abs_slip, max_dT, mean_rhoS_um2) on the
times common to all four runs:  e_k = ||h_k - h_ref||_L2(0,T) / ||h_ref||_L2(0,T)  for k = dt, dt/2, dt/4 (three data points), then
p_naive, p_corr (exact for steps halving each time, reference error included) and e_true_dt = e_1 + e_3/(2^p - 1), as in time_convergence_norms4.py.
A history that is identically zero (load_GPa and W_ext at 90 deg) gives nan.
Outputs in --out:  global_errors4.csv, global_errors4.json, tc_global4.png
"""
import argparse, json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import time_convergence as tc
import time_convergence_norms4 as n4


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--runs", nargs=4, required=True); ap.add_argument("--out", default="post/time_convergence_dt8")
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
    runs = sorted(a.runs, key=tc.dt_of, reverse=True); dts = [tc.dt_of(d) for d in runs]
    print("runs:", ", ".join("%s (dt=%g)" % (tc.tag_of(d), x) for d, x in zip(runs, dts)), " -> reference = last")
    if any(abs(dts[i] / dts[i + 1] - 2.0) > 1e-6 for i in range(3)): sys.exit("the four steps must halve each time: %s" % dts)
    L = [tc.read_log(d) for d in runs]
    tg = sorted(set.intersection(*[set(l) for l in L]))
    if len(tg) < 2: sys.exit("fewer than 2 common times in the logs")
    cols = [c for c in tc.GLOBAL if all(c in next(iter(l.values())) for l in L)]
    rows = []
    for c in cols:
        A = [np.array([float(l[t][c]) for t in tg]) for l in L]
        nrm = lambda x: float(np.sqrt(tc.trap(x ** 2, tg)))
        nr = nrm(A[3])
        rows.append(dict(quantity=c, norm_ref=nr, **n4.summarize(nrm(A[0] - A[3]), nrm(A[1] - A[3]), nrm(A[2] - A[3]), nr, nr),
                         end_dt=A[0][-1], end_dt2=A[1][-1], end_dt4=A[2][-1], end_dt8=A[3][-1]))
    tc.wcsv(os.path.join(a.out, "global_errors4.csv"), rows)
    json.dump(dict(runs=[tc.tag_of(d) for d in runs], dt=dts, n_times=len(tg), t_end=tg[-1], rows=rows),
              open(os.path.join(a.out, "global_errors4.json"), "w"), indent=1, default=float)
    print("global histories on %d common times, 0 .. %g ns; errors vs the dt/8 reference (dt = %g, dt/2, dt/4)" % (len(tg), tg[-1], dts[0]))
    print("  %-15s %9s %9s %9s  %7s %7s  %7s %7s  %10s" % ("quantity", "e(dt)", "e(dt/2)", "e(dt/4)", "p12 n", "p23 n", "p12 c", "p23 c", "e_true(dt)"))
    for r in rows:
        print("  %-15s %9.2e %9.2e %9.2e  %7.2f %7.2f  %7.2f %7.2f  %10.2e" % (r["quantity"], r["e_dt"], r["e_dt2"], r["e_dt4"],
              r["p_naive_12"], r["p_naive_23"], r["p_corr_12"], r["p_corr_23"], r["e_true_dt"]))
    ok = [r for r in rows if np.isfinite(r["e_true_dt"])]
    if ok:
        w = max(ok, key=lambda r: r["e_true_dt"])
        print("worst e_true(dt): %s %.2e  (%s 1 %%)" % (w["quantity"], w["e_true_dt"], "<=" if w["e_true_dt"] <= 0.01 else ">"))
    h = np.array(dts[:3]); fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    for i, r in enumerate(ok):
        ax[0].loglog(h, [r["e_dt"], r["e_dt2"], r["e_dt4"]], "o-", color="C%d" % (i % 10), lw=1, label=r["quantity"])
        ax[1].plot([1, 2], [r["p_corr_12"], r["p_corr_23"]], "o-", color="C%d" % (i % 10), lw=1)
    if ok:
        e0 = np.nanmedian([r["e_dt4"] for r in ok])
        for pp, ls in ((1, "k:"), (2, "k--")): ax[0].loglog(h, e0 * (h / h[-1]) ** pp, ls, lw=.8, label="slope %d" % pp)
    ax[0].set_xlabel("dt (ns)"); ax[0].set_ylabel("history error vs dt/8 (L2 in time, relative)"); ax[0].legend(fontsize=6, ncol=2); ax[0].grid(alpha=.3, which="both")
    ax[1].axhspan(0.6, 1.6, color="0.9"); ax[1].set_xticks([1, 2]); ax[1].set_xticklabels(["dt -> dt/2", "dt/2 -> dt/4"]); ax[1].set_ylabel("observed order (corrected)"); ax[1].grid(alpha=.3)
    fig.tight_layout(); fig.savefig(os.path.join(a.out, "tc_global4.png"), dpi=130); plt.close(fig)
    print("written to", a.out)


if __name__ == "__main__":
    main()
