#!/usr/bin/env python3
"""Fit the KWC grain-boundary energy, including the quadratic orientation term, to the Tschopp et al. (2015) symmetric tilt
energies (Cu <110> by default: the in-plane rotation axis of the model is [-110]).

Boundary energy of a 1D tilt boundary (kwc1d.py; regularisation delta, h(phi) = phi^2):
    gamma = f0 min int[(1-phi)^2/2 + a^2 phi'^2/2 + s g(phi)(sqrt(theta'^2+delta^2)-delta) + eps^2/2 h(phi) theta'^2] dx,
h(phi) = 1 (--h 1, default) or phi^2 (--h phi2)
gamma is linear in f0 (closed-form least squares); a sets the width and is held at --a; s and eps are fitted
(--eps X holds eps fixed).  The eps = 0 fit (sharp limit, fast) is reported alongside, with the production values.

usage:  python3 kwc_fit_tschopp.py tschopp_data/gb_energies.csv [--metal Cu] [--axis 110] [--max-angle 30]
            [--a 0.10] [--delta 1e-3] [--eps X | --eps-start 0.2] [--np N] [--out fit]
Runtime: about 1-3 s per boundary per evaluation and ~60-150 evaluations; use --np (one process per boundary)."""
import argparse, csv, math, os, sys
from multiprocessing import Pool
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kwc1d
from scipy.optimize import minimize_scalar, minimize

PROD = dict(f0=5.944e-3, a=0.10, s=1.05, eps=0.0)        # GPa, um, um, um
HFORM = "1"
PG = np.linspace(0, 1 - 1e-9, 200001); GG = -2*(PG + np.log1p(-PG))


def unit_gamma(args):
    """gamma for f0 = 1 GPa (gamma is proportional to f0)"""
    dth, a, s, eps, delta = args[:5]; hform = args[5] if len(args) > 5 else HFORM
    L = max(8*a, 8*eps, 0.8); n = int(max(1600, 2*L/0.001))
    return kwc1d.gamma(dth, 1.0, a, s, delta, eps=eps, L=L, n=n, hform=hform)[0]


def sharp_unit(th, a, s):
    return np.array([np.min(a*(1 - PG)**2 + s*t*GG) for t in th])*1e3


class Fitter:
    def __init__(self, th, E, a, delta, pool):
        self.th, self.E, self.a, self.delta, self.pool, self.n = th, E, a, delta, pool, 0
        self.hform = HFORM

    def unit(self, s, eps):
        self.n += 1
        return np.array(self.pool.map(unit_gamma, [(t, self.a, s, eps, self.delta, self.hform) for t in self.th]))

    def lsq(self, Gv):
        f0 = float(Gv @ self.E/(Gv @ Gv)); return f0, float(np.mean((f0*Gv - self.E)**2))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv"); ap.add_argument("--metal", default="Cu"); ap.add_argument("--axis", default="110")
    ap.add_argument("--max-angle", type=float, default=30.0); ap.add_argument("--a", type=float, default=0.10)
    ap.add_argument("--delta", type=float, default=1e-3); ap.add_argument("--eps", type=float, default=None)
    ap.add_argument("--eps-start", type=float, default=0.2)
    ap.add_argument("--np", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", "1")))
    ap.add_argument("--h", default="1", choices=("1", "phi2"), help="h(phi) of the eps^2/2 h |grad theta|^2 term")
    ap.add_argument("--out", default="fit")
    A = ap.parse_args()
    global HFORM; HFORM = A.h
    ang, en = [], []
    for r in csv.DictReader(open(A.csv)):
        if A.metal and r["metal"] and r["metal"].lower() != A.metal.lower(): continue
        if A.axis and r["axis"] and r["axis"] != A.axis: continue
        t = float(r["angle_deg"]) % 180.0; ang.append(min(t, 180.0 - t)); en.append(float(r["energy"]))
    ang, en = np.array(ang), np.array(en)
    if not len(ang): sys.exit("no rows for this metal/axis: check the metal/axis columns of the CSV")
    if np.median(en) > 10: en = en/1000.0; print("energies read as mJ/m^2 -> J/m^2")
    m = (ang > 0) & (ang <= A.max_angle); th = np.radians(ang[m]); E = en[m]
    print(f"{len(ang)} boundaries for {A.metal} <{A.axis}>; {m.sum()} with 0 < angle <= {A.max_angle} deg are fitted; {A.np} processes", flush=True)
    a = A.a; L = []
    # eps = 0, sharp limit (delta-independent to < 0.03 %)
    def r0(s):
        Gv = sharp_unit(th, a, s); f0 = float(Gv @ E/(Gv @ Gv)); return float(np.mean((f0*Gv - E)**2))
    s0 = minimize_scalar(lambda ls: r0(math.exp(ls)), bounds=(math.log(1e-3), math.log(1e2)), method="bounded").x; s0 = math.exp(s0)
    Gv = sharp_unit(th, a, s0); f00 = float(Gv @ E/(Gv @ Gv))
    L.append(f"eps = 0 fit:   f0 = {f00:.4e} GPa, s = {s0:.4f} um (a = {a} um);  rms {math.sqrt(r0(s0))*1e3:.1f} mJ/m^2")
    with Pool(A.np) as pool:
        F = Fitter(th, E, a, A.delta, pool)
        if A.eps is not None:
            def obj(ls):
                return F.lsq(F.unit(math.exp(ls), A.eps))[1]
            res = minimize_scalar(obj, bounds=(math.log(1e-3), math.log(1e2)), method="bounded", options=dict(xatol=1e-4))
            s1, e1 = math.exp(res.x), A.eps
        else:
            def obj(z):
                s, eps = math.exp(z[0]), math.exp(z[1]); v = F.lsq(F.unit(s, eps))[1]
                print(f"  eval {F.n}: s = {s:.4f}, eps = {eps:.4f} um, rms = {math.sqrt(v)*1e3:.2f} mJ/m^2", flush=True); return v
            res = minimize(obj, [math.log(max(s0, 1e-3)), math.log(A.eps_start)], method="Nelder-Mead",
                           options=dict(xatol=2e-3, fatol=1e-12, maxfev=200))
            s1, e1 = math.exp(res.x[0]), math.exp(res.x[1])
        Gv = F.unit(s1, e1); f01, v1 = F.lsq(Gv)
        L.append(f"eps fit:       f0 = {f01:.4e} GPa, s = {s1:.4f} um, eps = {e1:.4f} um (a = {a} um, delta = {A.delta}, h = {A.h});  rms {math.sqrt(v1)*1e3:.1f} mJ/m^2  [{F.n} evaluations]")
        Gp = F.unit(PROD["s"], PROD["eps"]); vp = float(np.mean((PROD["f0"]*Gp - E)**2))
        L.append(f"production:    f0 = {PROD['f0']:.4e} GPa, s = {PROD['s']} um, eps = 0 (a = {PROD['a']});  rms {math.sqrt(vp)*1e3:.1f} mJ/m^2")
        tt = np.radians([2, 5, 10, 15, 20, 30, 45])
        g_fit = f01*np.array(pool.map(unit_gamma, [(t, a, s1, e1, A.delta, A.h) for t in tt]))
        g_prod = PROD["f0"]*np.array(pool.map(unit_gamma, [(t, a, PROD["s"], 0.0, A.delta, A.h) for t in tt]))
    L.append("gamma (J/m^2)  angle: " + "  ".join(f"{d:5.0f}" for d in np.degrees(tt)))
    L.append("               fit:   " + "  ".join(f"{v:5.3f}" for v in g_fit))
    L.append("               prod:  " + "  ".join(f"{v:5.3f}" for v in g_prod))
    L.append("data used (angle deg, J/m^2): " + ", ".join(f"({x:.2f}, {y:.3f})" for x, y in zip(ang[m], E)))
    os.makedirs(A.out, exist_ok=True); open(os.path.join(A.out, "fit_summary.txt"), "w").write("\n".join(L) + "\n"); print("\n".join(L))
    try:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        plt.figure(figsize=(7, 4.5))
        plt.plot(ang, en, "o", color="0.65", ms=4, label="Tschopp et al. 2015 (all)")
        plt.plot(ang[m], E, "ko", ms=5, label="fitted")
        plt.plot(np.degrees(tt), g_fit, "C3o-", label=f"KWC fit: f0={f01:.3g}, s={s1:.3g}, eps={e1:.3g}")
        plt.plot(np.degrees(tt), g_prod, "C0s--", label="KWC production (eps = 0)")
        plt.xlabel("misorientation (deg)"); plt.ylabel("$\\gamma$ (J/m$^2$)"); plt.legend(fontsize=8); plt.grid(alpha=0.3)
        plt.tight_layout(); plt.savefig(os.path.join(A.out, "kwc_fit.png"), dpi=140)
    except ImportError:
        pass


if __name__ == "__main__":
    main()
