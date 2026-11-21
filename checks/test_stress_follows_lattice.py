#!/usr/bin/env python3
"""Stress follows the lattice: the elasticity tensor is rotated once (initial orientation) and never updated.

Reference: a crystal that starts at orientation alpha (C = R(alpha)*C^c), F = V, Fp = I.
Tests (C fixed at orientation 0; the lattice is turned by alpha during the history):
  T1 rigid rotation         F = V R(alpha),      Fp = I
  T2 grain-boundary spin    F = V,               Fp = R(-alpha)
  T3 slip + spin            Fp = slips on the three solver systems, then R(-alpha);  F = V R(alpha) Fp
  T4 time-stepped spin      200 steps of Fp <- exp(dw* What) Fp (as trialFp) to w* = -73 deg while F goes from I to V
  W  wrong on purpose       T2 with C rotated as well; must fail, and must equal a crystal at 2 alpha
Kirchhoff stress tau = P F^T, exact constraint; 35 stretch states x 12 orientations.
Pass: T1-T4 max rel. error < tol; W > 0.1 (the test discriminates) and W vs 2 alpha < tol.

usage: python3 checks/test_stress_follows_lattice.py [--law svk|prod|both] [--tol 1e-10] [--fig DIR | --no-fig]
exit 0 = pass, 1 = fail.  See checks/LATTICE_FRAME_TESTS.md.
"""
import argparse, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lattice_frame_common as L


def run(law, tol):
    V = L.load_set(); err = {k: 0.0 for k in ("T1", "T2", "T3", "T4", "W", "W2a")}; aniso = 0.0
    for Vk in V:
        t0 = L.state(Vk, np.eye(2), L.CC, law)["tau"]
        for al in L.ALPHAS:
            tref = L.state(Vk, np.eye(2), L.star(L.R3(al), L.CC), law)["tau"]; sc = np.abs(tref).max()
            aniso = max(aniso, np.abs(tref - t0).max()/sc)
            e = lambda t: np.abs(t - tref).max()/sc
            err["T1"] = max(err["T1"], e(L.state(Vk @ L.R2(al), np.eye(2), L.CC, law)["tau"]))
            err["T2"] = max(err["T2"], e(L.state(Vk, L.R2(-al), L.CC, law)["tau"]))
            Fp = L.history_Fp(al)
            err["T3"] = max(err["T3"], e(L.state(Vk @ L.R2(al) @ Fp, Fp, L.CC, law)["tau"]))
            tw = L.state(Vk, L.R2(-al), L.star(L.R3(al), L.CC), law)["tau"]
            err["W"] = max(err["W"], e(tw))
            t2 = L.state(Vk, np.eye(2), L.star(L.R3(2*al), L.CC), law)["tau"]
            err["W2a"] = max(err["W2a"], np.abs(tw - t2).max()/np.abs(t2).max())
    for Vk in V[::3]:
        for F, Fp, al in L.stepped_spin(Vk):
            tref = L.state(F, np.eye(2), L.star(L.R3(al), L.CC), law)["tau"]
            err["T4"] = max(err["T4"], np.abs(L.state(F, Fp, L.CC, law)["tau"] - tref).max()/np.abs(tref).max())
    ok = {k: err[k] < tol for k in ("T1", "T2", "T3", "T4", "W2a")}; ok["W"] = err["W"] > 0.1
    return err, ok, aniso, len(V)


def figure(law, path):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    Vf = np.diag([1.10, 0.95]); A = np.deg2rad(np.linspace(-90, 90, 361)); As = np.deg2rad(np.arange(-90, 91, 10))
    tau = lambda F, Fp, C: L.state(F, Fp, C, law)["tau"]
    ref = np.array([tau(Vf, np.eye(2), L.star(L.R3(a), L.CC)) for a in A])
    tests = [("T1: rigid rotation, $\\mathbb{C}$ fixed\n$F=V R(\\alpha)$, $F_p=I$", lambda a: tau(Vf @ L.R2(a), np.eye(2), L.CC)),
             ("T2: grain-boundary spin, $\\mathbb{C}$ fixed\n$F=V$, $F_p=R(-\\alpha)$", lambda a: tau(Vf, L.R2(-a), L.CC)),
             ("T3: three slips + spin in $F_p$, $\\mathbb{C}$ fixed\n$F=V R(\\alpha) F_p$", lambda a: tau(Vf @ L.R2(a) @ L.history_Fp(a), L.history_Fp(a), L.CC)),
             None,
             ("W: spin AND $\\mathbb{C}$ rotated (wrong)\n$F=V$, $F_p=R(-\\alpha)$, $\\mathbb{C}=R(\\alpha)\\star\\mathbb{C}^c$", lambda a: tau(Vf, L.R2(-a), L.star(L.R3(a), L.CC)))]
    comp = [((0, 0), "$\\tau_{11}$", "C0"), ((1, 1), "$\\tau_{22}$", "C1"), ((0, 1), "$\\tau_{12}$", "C3")]
    fig, axs = plt.subplots(2, 3, figsize=(15, 9)); axs = axs.ravel()
    for k, t in enumerate(tests):
        ax = axs[k]
        if t is None:
            ts, tt, tr = [], [], []
            for n, (F, Fp, al) in enumerate(L.stepped_spin(Vf, every=5)):
                ts.append(5*(n + 1)); tt.append(tau(F, Fp, L.CC)); tr.append(tau(F, np.eye(2), L.star(L.R3(al), L.CC)))
            ts, tt, tr = np.array(ts), np.array(tt), np.array(tr)
            for (i, j), nm, c in comp:
                ax.plot(ts, tr[:, i, j], "-", color=c, lw=1.6, label=nm + " reference"); ax.plot(ts[1::2], tt[1::2, i, j], "o", color=c, mfc="none", ms=6)
            ax.set_xlabel("step (lattice angle grows linearly to 73 deg; F: I to V)")
            ax.set_title("T4: grain-boundary spin stepped in time,\n$\\mathbb{C}$ fixed, 200 steps", fontsize=10.5)
            ax.legend(fontsize=8.5, loc="upper left"); ev = np.abs(tt - tr).max()/np.abs(tr).max()
        else:
            tt = np.array([t[1](a) for a in As]); tra = np.array([tau(Vf, np.eye(2), L.star(L.R3(a), L.CC)) for a in As])
            for (i, j), nm, c in comp:
                ax.plot(np.rad2deg(A), ref[:, i, j], "-", color=c, lw=1.6, label=nm + " reference")
                ax.plot(np.rad2deg(As), tt[:, i, j], "x" if k == 4 else "o", color=c, mfc="none", ms=6)
            if k == 4:
                r2 = np.array([tau(Vf, np.eye(2), L.star(L.R3(2*a), L.CC)) for a in A])
                ax.plot(np.rad2deg(A), r2[:, 0, 0], ":", color="C0", lw=1.2, label="$\\tau_{11}$ of a crystal at $2\\alpha$"); ax.legend(fontsize=8, loc="upper right")
            ax.set_xlabel("lattice rotation $\\alpha$ (deg)"); ax.set_xticks(range(-90, 91, 30)); ax.set_title(t[0], fontsize=10.5)
            ev = np.abs(tt - tra).max()/np.abs(ref).max()
        ax.set_ylabel("Kirchhoff stress $\\tau=PF^T$ (GPa)")
        ax.text(0.02, 0.03, f"max rel. error vs reference: {ev:.1e}", transform=ax.transAxes, fontsize=9.5,
                color="firebrick" if k == 4 else "darkgreen", bbox=dict(fc="white", ec="none", alpha=0.8))
    return fig, axs[5]


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--law", default="both", choices=("svk", "prod", "both")); ap.add_argument("--tol", type=float, default=1e-10)
    ap.add_argument("--fig", default="checks/out", help="folder for the figure (default checks/out)"); ap.add_argument("--no-fig", action="store_true")
    a = ap.parse_args()
    src, msgs = L.check_against_solver(); bad = False
    print("solver constants:", src or "-")
    for lvl, m in msgs:
        print(f"  [{lvl}] {m}"); bad |= lvl == "FAIL"
    laws = ("svk", "prod") if a.law == "both" else (a.law,)
    for law in laws:
        err, ok, aniso, nV = run(law, a.tol)
        print(f"\nstress follows the lattice, law = {law}: {nV} stretch states x {len(L.ALPHAS)} orientations; "
              f"orientation changes the stress by up to {100*aniso:.0f} %")
        lab = {"T1": "T1 rigid rotation, C fixed", "T2": "T2 grain-boundary spin, C fixed", "T3": "T3 slip + spin, C fixed",
               "T4": "T4 time-stepped spin, C fixed", "W": "W  spin AND C rotated (must differ)", "W2a": "W  equals a crystal at 2 alpha"}
        for k in ("T1", "T2", "T3", "T4", "W", "W2a"):
            print(f"  {lab[k]:38s} {err[k]:9.2e}  {'PASS' if ok[k] else 'FAIL'}"); bad |= not ok[k]
        if not a.no_fig:
            try:
                fig, ax = figure(law, None)
                keys = ["T1", "T2", "T3", "T4", "W"]
                ax.bar(keys, [err[k] for k in keys], color=["C2"]*4 + ["C3"]); ax.set_yscale("log"); ax.set_ylim(1e-16, 10)
                for i, k in enumerate(keys): ax.text(i, err[k]*2.2, f"{err[k]:.1e}", ha="center", fontsize=9)
                ax.axhline(a.tol, color="k", ls="--", lw=0.8); ax.text(4.45, a.tol*1.6, "tolerance", ha="right", fontsize=8.5)
                ax.set_ylabel(f"max rel. error, {nV} stretches x {len(L.ALPHAS)} orientations"); ax.set_title("Summary over the full test set", fontsize=10.5)
                fig.suptitle(f"Stress follows the lattice ({law} law). Lines: crystal starting at $\\alpha$. Markers: tests. "
                             "$V=\\mathrm{diag}(1.10,0.95)$, Cu, plane strain, $\\mu_c=754$ GPa, exact constraint", fontsize=11)
                fig.tight_layout(rect=(0, 0, 1, 0.96)); os.makedirs(a.fig, exist_ok=True)
                p = os.path.join(a.fig, f"stress_follows_lattice_{law}.png"); fig.savefig(p, dpi=140); print("  figure:", p)
            except ImportError:
                print("  (matplotlib not available: no figure)")
    print("\nRESULT:", "FAIL" if bad else "PASS"); sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
