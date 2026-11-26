#!/usr/bin/env python3
"""time_convergence.py  --  time-step convergence of the production set-up from three runs at dt, dt/2, dt/4.

usage:  python3 time_convergence.py --runs out_<tag dt> out_<tag dt/2> out_<tag dt/4> [--out post/time_convergence]
                                    [--every 0.5] [--grain-times 3 6 9 12 15] [--np 4]

The three runs share the mesh (Lagrangian, same node and element order) and the relaxed state, so fields are compared node by
node / element by element at equal times (times from <tag>.pvd; frames keep their spacing in ns at every dt).

Error measures (q_1, q_2, q_4 = solutions at dt, dt/2, dt/4; the dt/4 run is the reference):
  global histories (from <tag>_log.csv, on the dt time grid):
      e_k = ||a_k - a_4||_L2(0,T) / ||a_4||_L2(0,T)                     (k = 1, 2)
  fields at each compared time t (area-weighted L2 over the domain; nodal fields with lumped P1 weights, element fields with
  the triangle areas), relative to the change since the relaxed state, so a field that barely changed is not called converged
  because of a large offset:
      e_k(t) = ||f_k(t) - f_4(t)|| / ||f_4(t) - f_4(0)||
  observed order (needs all three):          p = log2( ||q_1 - q_2|| / ||q_2 - q_4|| )
  Richardson error of the production run:    E_1 = ||q_1 - q_2|| / (1 - 2^-p)   (p clipped to [0.5, 2]; also given with p = 1)
                                             relative to the same norm as e_k.

Acceptance criteria for the production dt (summary.json, PASS/FAIL each):
  C1  order:   observed p of the global histories in [0.6, 1.6] (backward Euler in the local, phase-field and heat updates is first
               order, so p ~ 1 is the asymptotic regime; outside it the Richardson estimate is not trusted).
  C2  global:  Richardson error E_1 <= 1 % for every global history (load, energies, dissipation, heat, mean slip, max dT).
  C3  fields:  Richardson error E_1 <= 5 % at every compared time for theta_lattice, phi, T, u, v and the accumulated slip
               (median over time also reported; the worst time is usually the impact or a localization event).
  C4  grains:  at every grain time, the grain count of the dt run within max(1, 10 %) of the dt/4 run, and its low- and high-angle
               boundary length fractions within 0.02 (absolute) of the dt/4 run (grain_count.py, GRAD/MERGE/MINAREA as grain_post,
               SYM 180).
  C3 failing while C1, C2 and C4 pass means the pointwise fields are path-sensitive (localization, boundary nucleation sites)
  while the integrated and statistical response is converged; summary.json says so.

Outputs in --out:  global_errors.csv, field_errors.csv, grain_stats.csv, summary.json, tc_global.png, tc_fields.png, tc_grains.png
"""
import argparse, csv, json, os, re, sys
from multiprocessing import Pool
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import grain_count as gc

E = lambda k, d: os.environ.get(k, d)
GLOBAL = ["load_GPa", "KE_trans", "E_elastic", "E_stored", "Q_heat", "W_ext", "D_cum", "mean_abs_slip", "max_dT",
          "mean_rhoS_um2"]
NODAL = ["theta_lattice_deg", "phi", "T", "u", "v"]   # u displacement, v velocity (vectors)
CELL = ["slip_path_phys", "rhoS", "rhoG", "wstar", "tau_c"]
FIELDS_C3 = ["theta_lattice_deg", "phi", "T", "u", "v", "slip_path_phys"]


def tag_of(d):
    return re.sub(r"^out_", "", os.path.basename(os.path.normpath(d)))


def dt_of(d):
    a = open(os.path.join(d, tag_of(d) + ".args")).read()
    return float(re.search(r"-dt ([^ ]+)", a).group(1))


def read_log(d):
    f = os.path.join(d, tag_of(d) + "_log.csv")
    rows = list(csv.DictReader(open(f)))
    out = {}
    for r in rows:
        try:
            t = round(float(r["t_ns"]), 6)
        except (KeyError, ValueError):
            continue
        out[t] = r                                               # later rows (after a resume) replace earlier ones
    return out


def frames(d):
    b = tag_of(d); out = {}
    rel = os.path.join(d, b + "_relaxed.vtu")
    if os.path.exists(rel):
        out[0.0] = rel
    pvd = os.path.join(d, b + ".pvd")
    for t, f in re.findall(r'timestep="([^"]+)" file="([^"]+)"', open(pvd).read()):
        p = os.path.join(d, os.path.basename(f))
        if os.path.exists(p) and "relaxed" not in p:
            out[round(float(t), 6)] = p
    return out


def load_fields(fn):
    v = gc.read_vtu(fn)
    pts = v["points"][:, :2].astype(float)
    conn, off, typ = v["connectivity"], v["offsets"], v["types"]
    start = np.r_[0, off[:-1]]
    tri = np.array([conn[s:s + 3] for s, t in zip(start, typ) if t == 5])
    istri = typ == 5
    F = {}
    for k in NODAL:
        a = v[k].astype(float)
        F[k] = a[:, :2] if a.ndim == 2 else a.ravel()             # u (displacement) and v (velocity) are vectors
    for k in CELL:
        F[k] = v[k].astype(float).ravel()[istri]
    return pts, tri, F


def weights(pts, tri):
    p0, p1, p2 = pts[tri[:, 0]], pts[tri[:, 1]], pts[tri[:, 2]]
    A = 0.5 * np.abs((p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1]) - (p2[:, 0] - p0[:, 0]) * (p1[:, 1] - p0[:, 1]))
    wn = np.bincount(tri.ravel(), weights=np.repeat(A / 3, 3), minlength=len(pts))
    return wn, A


def trap(y, x):
    x = np.asarray(x, float); y = np.asarray(y, float)
    return float(np.sum(0.5 * (y[1:] + y[:-1]) * np.diff(x)))


def l2(x, w):
    return float(np.sqrt(np.sum((w[:, None] if x.ndim == 2 else w) * x * x)))


def order(d12, d24):
    if d12 <= 0 or d24 <= 0:
        return float("nan")
    return float(np.log2(d12 / d24))


def richardson(d12, p):
    pc = 1.0 if not np.isfinite(p) else min(max(p, 0.5), 2.0)
    return d12 / (1 - 2.0 ** -pc), d12 / 0.5


def field_job(args):
    t, f1, f2, f4, f40 = args
    P, T, F1 = load_fields(f1); _, _, F2 = load_fields(f2); _, _, F4 = load_fields(f4); _, _, F0 = load_fields(f40)
    wn, wc = weights(P, T)
    rows = []
    for k in NODAL + CELL:
        w = wn if k in NODAL else wc
        ref = l2(F4[k] - F0[k], w)
        d14, d24, d12 = l2(F1[k] - F4[k], w), l2(F2[k] - F4[k], w), l2(F1[k] - F2[k], w)
        p = order(d12, d24)
        R, R1 = richardson(d12, p)
        s = ref if ref > 0 else float("nan")
        rows.append(dict(t_ns=t, field=k, ref_change=ref, e_dt=d14 / s, e_dt2=d24 / s, p=p, rich_dt=R / s, rich_dt_p1=R1 / s))
    return rows


def grain_job(args):
    t, fs, par = args
    out = []
    for lab, f in fs:
        pts, tri, th, phi, label = gc.load(f)
        M = gc.Mesh(pts, tri, th, phi)
        seg = gc.segment(M, par["grad"], par["minarea"], par["merge"])
        if seg is None:
            out.append(dict(t_ns=t, run=lab, grains=0, L=0, f_lagb=0, f_hagb=0)); continue
        rows, bnd = gc.analyse(M, seg, label, 180.0)
        L = sum(b["length_um"] for b in bnd)
        lo = sum(b["length_um"] for b in bnd if 2 <= b["misorient_sym_deg"] < 15)
        hi = sum(b["length_um"] for b in bnd if b["misorient_sym_deg"] >= 15)
        out.append(dict(t_ns=t, run=lab, grains=seg["G"], L=L, f_lagb=lo / L if L else 0, f_hagb=hi / L if L else 0))
    return out


def wcsv(fn, rows):
    with open(fn, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--runs", nargs=3, required=True, help="run folders at dt, dt/2, dt/4 (any order: sorted by dt)")
    ap.add_argument("--out", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "time_convergence"))
    ap.add_argument("--every", type=float, default=0.5, help="compare fields every this many ns (0.5)")
    ap.add_argument("--grain-times", type=float, nargs="*", default=[3, 6, 9, 12, 15])
    ap.add_argument("--np", type=int, default=int(E("NP", 4)))
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    runs = sorted(a.runs, key=dt_of, reverse=True)                 # dt, dt/2, dt/4
    dts = [dt_of(d) for d in runs]
    labs = [f"dt={x:g}" for x in dts]
    print("runs:", ", ".join(f"{tag_of(d)} ({l})" for d, l in zip(runs, labs)))
    for i in (1, 2):
        if abs(dts[0] / dts[i] - 2 ** i) > 1e-6:
            print(f"WARNING: dt ratios are not 1 : 1/2 : 1/4 ({dts}); the order formula assumes them")
    S = dict(runs=[tag_of(d) for d in runs], dt=dts, criteria={})

    # ---- global histories
    L = [read_log(d) for d in runs]
    tg = sorted(set(L[0]) & set(L[1]) & set(L[2]))
    S["t_end_common"] = tg[-1] if tg else None
    cols = [c for c in GLOBAL if all(c in next(iter(l.values())) for l in L)]
    grows, H = [], {}
    for c in cols:
        A = [np.array([float(l[t][c]) for t in tg]) for l in L]
        H[c] = A
        n4 = np.sqrt(trap(A[2] ** 2, tg))
        nrm = lambda x: np.sqrt(trap(x ** 2, tg))
        d14, d24, d12 = nrm(A[0] - A[2]), nrm(A[1] - A[2]), nrm(A[0] - A[1])
        p = order(d12, d24); R, R1 = richardson(d12, p)
        s = n4 if n4 > 0 else float("nan")
        grows.append(dict(quantity=c, norm_ref=n4, e_dt=d14 / s, e_dt2=d24 / s, p=p, rich_dt=R / s, rich_dt_p1=R1 / s,
                          end_dt=A[0][-1], end_dt2=A[1][-1], end_dt4=A[2][-1]))
    wcsv(os.path.join(a.out, "global_errors.csv"), grows)
    print(f"global histories on {len(tg)} common times, 0 .. {S['t_end_common']} ns")
    for r in grows:
        print(f"  {r['quantity']:<15} e(dt) {r['e_dt']:.2e}  e(dt/2) {r['e_dt2']:.2e}  p {r['p']:5.2f}  Richardson E(dt) {r['rich_dt']:.2e}")

    # ---- fields
    Fr = [frames(d) for d in runs]
    tf = sorted(set(Fr[0]) & set(Fr[1]) & set(Fr[2]))
    have0 = all(0.0 in f for f in Fr)
    if have0:
        t0 = 0.0
    else:
        c = [t for t in tf if t > 0]
        if len(c) < 2:
            sys.exit("no relaxed frame and fewer than 2 common frames; frames per run: %s" % [len(f) for f in Fr])
        t0 = c[0]
        print("note: no relaxed frame in the run folders; reference state = the finest-dt frame at t = %g ns (denominators: change since then; errors start after it)" % t0)
    pick = [t for t in tf if t > t0 and abs(t / a.every - round(t / a.every)) < 1e-6]
    if tf and tf[-1] not in pick and tf[-1] > t0:
        pick.append(tf[-1])
    jobs = [(t, Fr[0][t], Fr[1][t], Fr[2][t], Fr[2][t0]) for t in pick]
    with Pool(a.np) as pool:
        frows = [r for rr in pool.map(field_job, jobs) for r in rr]
    wcsv(os.path.join(a.out, "field_errors.csv"), frows)
    print(f"fields at {len(pick)} times")

    # ---- grains
    par = dict(grad=float(E("GRAD", 2)), merge=float(E("MERGE", 2)), minarea=float(E("MINAREA", 0.1)))
    gt = [t for t in a.grain_times if all(round(t, 6) in f for f in Fr)]
    with Pool(min(a.np, max(len(gt), 1))) as pool:
        grows2 = [r for rr in pool.map(grain_job, [(t, [(l, f[round(t, 6)]) for l, f in zip(labs, Fr)], par) for t in gt]) for r in rr]
    if grows2:
        wcsv(os.path.join(a.out, "grain_stats.csv"), grows2)

    # ---- criteria
    ps = [r["p"] for r in grows if np.isfinite(r["p"])]
    S["criteria"]["C1_order"] = dict(p_global=dict((r["quantity"], r["p"]) for r in grows),
                                     pass_=bool(ps) and all(0.6 <= p <= 1.6 for p in ps))
    S["criteria"]["C2_global"] = dict(worst=max(grows, key=lambda r: r["rich_dt"])["quantity"] if grows else None,
                                      worst_rich=max((r["rich_dt"] for r in grows), default=float("nan")),
                                      pass_=all(r["rich_dt"] <= 0.01 for r in grows))
    c3 = {}
    for k in FIELDS_C3:
        v = np.array([r["rich_dt"] for r in frows if r["field"] == k])
        v = v[np.isfinite(v)]
        c3[k] = dict(max=float(v.max()) if len(v) else None, median=float(np.median(v)) if len(v) else None,
                     t_worst=[r["t_ns"] for r in frows if r["field"] == k and r["rich_dt"] == v.max()][0] if len(v) else None)
    S["criteria"]["C3_fields"] = dict(fields=c3, pass_=all(c["max"] is not None and c["max"] <= 0.05 for c in c3.values()))
    ok4, det = True, []
    for t in gt:
        g = {r["run"]: r for r in grows2 if r["t_ns"] == t}
        a1, a4 = g[labs[0]], g[labs[2]]
        okc = abs(a1["grains"] - a4["grains"]) <= max(1, 0.1 * a4["grains"])
        okf = abs(a1["f_lagb"] - a4["f_lagb"]) <= 0.02 and abs(a1["f_hagb"] - a4["f_hagb"]) <= 0.02
        ok4 &= okc and okf
        det.append(dict(t_ns=t, grains=[g[l]["grains"] for l in labs], f_lagb=[g[l]["f_lagb"] for l in labs],
                        f_hagb=[g[l]["f_hagb"] for l in labs], pass_=okc and okf))
    S["criteria"]["C4_grains"] = dict(times=det, pass_=ok4 and bool(det))
    C = S["criteria"]
    S["verdict"] = ("production dt converged" if all(C[k]["pass_"] for k in C) else
                    "integrated and statistical response converged, pointwise fields path-sensitive"
                    if C["C1_order"]["pass_"] and C["C2_global"]["pass_"] and C["C4_grains"]["pass_"] else
                    "production dt NOT converged by these criteria")
    json.dump(S, open(os.path.join(a.out, "summary.json"), "w"), indent=1, default=float)
    for k, c in C.items():
        print(f"{k}: {'PASS' if c['pass_'] else 'FAIL'}")
    print("verdict:", S["verdict"])

    # ---- figures
    sel = [c for c in ["load_GPa", "KE_trans", "E_elastic", "D_cum", "Q_heat", "max_dT"] if c in H]
    if sel:
        fig, ax = plt.subplots(2, len(sel), figsize=(3.2 * len(sel), 5.5), sharex=True, squeeze=False)
        for j, c in enumerate(sel):
            for i, (A, l) in enumerate(zip(H[c], labs)):
                ax[0, j].plot(tg, A, lw=1, label=l, color=f"C{i}")
            ax[1, j].plot(tg, H[c][0] - H[c][2], lw=1, color="C0", label=f"{labs[0]} − {labs[2]}")
            ax[1, j].plot(tg, H[c][1] - H[c][2], lw=1, color="C1", label=f"{labs[1]} − {labs[2]}")
            ax[0, j].set_title(c, fontsize=9); ax[1, j].set_xlabel("t (ns)")
        ax[0, 0].legend(fontsize=7); ax[1, 0].legend(fontsize=7); ax[1, 0].set_ylabel("difference")
        fig.tight_layout(); fig.savefig(os.path.join(a.out, "tc_global.png"), dpi=130); plt.close(fig)
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    for i, k in enumerate(NODAL + CELL):
        rr = [r for r in frows if r["field"] == k]
        t = [r["t_ns"] for r in rr]
        ax[0].semilogy(t, [r["e_dt"] for r in rr], "-", color=f"C{i}", lw=1, label=k)
        ax[0].semilogy(t, [r["e_dt2"] for r in rr], "--", color=f"C{i}", lw=1)
        ax[1].plot(t, [r["p"] for r in rr], "-", color=f"C{i}", lw=1, label=k)
    ax[0].axhline(0.05, color="k", lw=.6, ls=":"); ax[0].set_title(f"field error vs {labs[2]}: solid {labs[0]}, dashed {labs[1]}", fontsize=9)
    ax[1].axhspan(0.6, 1.6, color="0.9"); ax[1].set_ylim(-0.5, 3); ax[1].set_title("observed order p", fontsize=9)
    for x in ax: x.set_xlabel("t (ns)"); x.grid(alpha=.3)
    ax[0].legend(fontsize=6, ncol=2); fig.tight_layout(); fig.savefig(os.path.join(a.out, "tc_fields.png"), dpi=130); plt.close(fig)
    if grows2:
        fig, ax = plt.subplots(1, 3, figsize=(11, 3.4))
        for i, l in enumerate(labs):
            rr = [r for r in grows2 if r["run"] == l]
            for j, k in enumerate(["grains", "f_lagb", "f_hagb"]):
                ax[j].plot([r["t_ns"] for r in rr], [r[k] for r in rr], "o-", color=f"C{i}", label=l)
        for x, k in zip(ax, ["grains", "LAGB length fraction (2–15°)", "HAGB length fraction (≥15°)"]):
            x.set_title(k, fontsize=9); x.set_xlabel("t (ns)"); x.grid(alpha=.3)
        ax[0].legend(fontsize=7); fig.tight_layout(); fig.savefig(os.path.join(a.out, "tc_grains.png"), dpi=130); plt.close(fig)
    print("written to", a.out)


if __name__ == "__main__":
    main()
