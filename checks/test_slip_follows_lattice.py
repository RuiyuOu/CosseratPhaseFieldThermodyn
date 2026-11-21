#!/usr/bin/env python3
"""Slip follows the lattice: the slip directions s_a and normals n_a are rotated once (initial orientation) and never updated.

Driving force: tau_a = M : (s_a x n_a) = tau : (Fe s_a x Fe^-T n_a), the stress resolved on the CURRENT lattice slip system.
Kinematics:    Fe Lp Fe^-1 = sum gdot_a (Fe s_a) x (Fe^-T n_a), slip along the current direction on the current plane.
Reference: a crystal that starts at orientation alpha (C and the three solver slip systems rotated by alpha), F = V, Fp = I.
Tests (C and the slip systems fixed at orientation 0; the lattice is turned by alpha during the history):
  T1 rigid rotation, T2 grain-boundary spin, T3 slip + spin, T4 time-stepped spin (as in test_stress_follows_lattice.py)
  W  wrong on purpose: T2 with C and the slip systems rotated as well (must fail; equals a crystal at 2 alpha)
  Ws wrong on purpose: T2 with only the slip systems rotated (must fail)
Compared: tau_a for the three systems, the current slip directions Fe s_a and normals Fe^-T n_a, and the push-forward
identity of the plastic velocity gradient.  35 stretch states x 12 orientations.

usage: python3 checks/test_slip_follows_lattice.py [--law svk|prod|both] [--tol 1e-10] [--fig DIR | --no-fig]
exit 0 = pass, 1 = fail.  See checks/LATTICE_FRAME_TESTS.md.
"""
import argparse, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lattice_frame_common as L

KEYS = ("T1", "T2", "T3", "T4", "W", "Ws")


def compare(test, ref, sc):
    et = max(abs(a[0] - b[0]) for a, b in zip(test, ref))/sc
    es = max(np.abs(a[1] - b[1]).max() for a, b in zip(test, ref))
    en = max(np.abs(a[2] - b[2]).max() for a, b in zip(test, ref))
    return et, es, en


def run(law, tol):
    V = L.load_set(); S0 = L.slip_systems(0.0)
    et = {k: 0.0 for k in KEYS}; ev = {k: 0.0 for k in KEYS}; e2a = 0.0; epf = 0.0
    def upd(k, r):
        et[k] = max(et[k], r[0]); ev[k] = max(ev[k], r[1], r[2])
    for Vk in V:
        for al in L.ALPHAS:
            Sa = L.slip_systems(al)
            st = L.state(Vk, np.eye(2), L.star(L.R3(al), L.CC), law); ref = L.resolved(st, Sa); sc = np.abs(st["tau"]).max()
            upd("T1", compare(L.resolved(L.state(Vk @ L.R2(al), np.eye(2), L.CC, law), S0), ref, sc))
            upd("T2", compare(L.resolved(L.state(Vk, L.R2(-al), L.CC, law), S0), ref, sc))
            Fp = L.history_Fp(al)
            upd("T3", compare(L.resolved(L.state(Vk @ L.R2(al) @ Fp, Fp, L.CC, law), S0), ref, sc))
            w = L.resolved(L.state(Vk, L.R2(-al), L.star(L.R3(al), L.CC), law), Sa)
            upd("W", compare(w, ref, sc))
            st2 = L.state(Vk, np.eye(2), L.star(L.R3(2*al), L.CC), law)
            e2a = max(e2a, compare(w, L.resolved(st2, L.slip_systems(2*al)), np.abs(st2["tau"]).max())[0])
            upd("Ws", compare(L.resolved(L.state(Vk, L.R2(-al), L.CC, law), Sa), ref, sc))
            # push-forward of the plastic velocity gradient onto the current lattice
            Fe = L.state(Vk @ L.R2(al) @ Fp, Fp, L.CC, law)["Fe"]; g = (0.3, -0.7, 1.1)
            Lp = sum(gi*np.outer(s, n) for gi, (s, n) in zip(g, S0))
            rhs = sum(gi*np.outer(Fe @ s, np.linalg.inv(Fe).T @ n) for gi, (s, n) in zip(g, S0))
            epf = max(epf, np.abs(Fe @ Lp @ np.linalg.inv(Fe) - rhs).max())
    for Vk in V[::3]:
        for F, Fp, al in L.stepped_spin(Vk):
            st = L.state(F, np.eye(2), L.star(L.R3(al), L.CC), law)
            upd("T4", compare(L.resolved(L.state(F, Fp, L.CC, law), S0), L.resolved(st, L.slip_systems(al)), np.abs(st["tau"]).max()))
    ok = {k: (et[k] < tol and ev[k] < tol) for k in ("T1", "T2", "T3", "T4")}
    ok["W"] = et["W"] > 0.05; ok["Ws"] = et["Ws"] > 0.05; ok["W2a"] = e2a < tol; ok["pf"] = epf < tol
    return et, ev, e2a, epf, ok, len(V)


def figure(law):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    Vf = np.diag([1.10, 0.95]); A = np.deg2rad(np.linspace(-90, 90, 361)); As = np.deg2rad(np.arange(-90, 91, 10)); S0 = L.slip_systems(0.0)
    taus = lambda F, Fp, C, S: np.array([r[0] for r in L.resolved(L.state(F, Fp, C, law), S)])
    refs = lambda a: taus(Vf, np.eye(2), L.star(L.R3(a), L.CC), L.slip_systems(a))
    ref = np.array([refs(a) for a in A])
    tests = [("T1: rigid rotation, $\\mathbb{C}$, $s_\\alpha$, $n_\\alpha$ fixed\n$F=V R(\\alpha)$, $F_p=I$", lambda a: taus(Vf @ L.R2(a), np.eye(2), L.CC, S0)),
             ("T2: grain-boundary spin, slip systems fixed\n$F=V$, $F_p=R(-\\alpha)$", lambda a: taus(Vf, L.R2(-a), L.CC, S0)),
             ("T3: three slips + spin in $F_p$, slip systems fixed\n$F=V R(\\alpha) F_p$", lambda a: taus(Vf @ L.R2(a) @ L.history_Fp(a), L.history_Fp(a), L.CC, S0)),
             None,
             ("W: spin AND $\\mathbb{C}$, slip systems rotated (wrong)", lambda a: taus(Vf, L.R2(-a), L.star(L.R3(a), L.CC), L.slip_systems(a))),
             ("Ws: spin AND slip systems rotated, $\\mathbb{C}$ fixed (wrong)", lambda a: taus(Vf, L.R2(-a), L.CC, L.slip_systems(a)))]
    cols = ["C0", "C1", "C3"]; names = ["system 1 $(1,0)$", "system 2 $(\\cos a_R, \\sin a_R)$", "system 3 $(\\cos a_R, -\\sin a_R)$"]
    fig, axs = plt.subplots(2, 4, figsize=(20, 9)); axs = axs.ravel()
    for k, t in enumerate(tests):
        ax = axs[k]; wrong = k >= 4
        if t is None:
            ts, tt, tr = [], [], []
            for n, (F, Fp, al) in enumerate(L.stepped_spin(Vf, every=5)):
                ts.append(5*(n + 1)); tt.append(taus(F, Fp, L.CC, S0)); tr.append(taus(F, np.eye(2), L.star(L.R3(al), L.CC), L.slip_systems(al)))
            ts, tt, tr = np.array(ts), np.array(tt), np.array(tr)
            for i in range(3):
                ax.plot(ts, tr[:, i], "-", color=cols[i], lw=1.6, label=names[i] + " reference"); ax.plot(ts[1::2], tt[1::2, i], "o", color=cols[i], mfc="none", ms=6)
            ax.set_xlabel("step (lattice angle grows linearly to 73 deg; F: I to V)"); ax.legend(fontsize=8, loc="upper left")
            ax.set_title("T4: grain-boundary spin stepped in time,\nslip systems fixed, 200 steps", fontsize=10.5)
            ev = np.abs(tt - tr).max()/np.abs(tr).max()
        else:
            tt = np.array([t[1](a) for a in As]); tra = np.array([refs(a) for a in As])
            for i in range(3):
                ax.plot(np.rad2deg(A), ref[:, i], "-", color=cols[i], lw=1.6, label=names[i] + " reference")
                ax.plot(np.rad2deg(As), tt[:, i], "x" if wrong else "o", color=cols[i], mfc="none", ms=6)
            ax.set_xlabel("lattice rotation $\\alpha$ (deg)"); ax.set_xticks(range(-90, 91, 30)); ax.set_title(t[0], fontsize=10.5)
            ev = np.abs(tt - tra).max()/np.abs(ref).max()
        ax.set_ylabel("resolved shear stress $\\tau_\\alpha = M:(s_\\alpha\\otimes n_\\alpha)$ (GPa)")
        ax.text(0.02, 0.03, f"max rel. error vs reference: {ev:.1e}", transform=ax.transAxes, fontsize=9.5,
                color="firebrick" if wrong else "darkgreen", bbox=dict(fc="white", ec="none", alpha=0.8))
    return fig, axs[6], axs[7]


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
        et, ev, e2a, epf, ok, nV = run(law, a.tol)
        print(f"\nslip follows the lattice, law = {law}: {nV} stretch states x {len(L.ALPHAS)} orientations, 3 slip systems")
        print(f"  {'case':44s} {'tau_a':>9s} {'Fe s, Fe^-T n':>14s}")
        lab = {"T1": "T1 rigid rotation, slip systems fixed", "T2": "T2 grain-boundary spin, slip systems fixed",
               "T3": "T3 slip + spin, slip systems fixed", "T4": "T4 time-stepped spin, slip systems fixed",
               "W": "W  spin AND C, s, n rotated (must differ)", "Ws": "Ws spin AND s, n rotated (must differ)"}
        for k in KEYS:
            print(f"  {lab[k]:44s} {et[k]:9.2e} {ev[k]:14.2e}  {'PASS' if ok[k] else 'FAIL'}"); bad |= not ok[k]
        print(f"  {'W  equals a crystal at 2 alpha':44s} {e2a:9.2e} {'':14s}  {'PASS' if ok['W2a'] else 'FAIL'}"); bad |= not ok["W2a"]
        print(f"  {'Fe Lp Fe^-1 = sum gdot (Fe s) x (Fe^-T n)':44s} {epf:9.2e} {'':14s}  {'PASS' if ok['pf'] else 'FAIL'}"); bad |= not ok["pf"]
        if not a.no_fig:
            try:
                fig, ax1, ax2 = figure(law)
                ax1.bar(list(KEYS), [et[k] for k in KEYS], color=["C2"]*4 + ["C3"]*2); ax1.set_yscale("log"); ax1.set_ylim(1e-16, 10)
                for i, k in enumerate(KEYS): ax1.text(i, et[k]*2.2, f"{et[k]:.1e}", ha="center", fontsize=8.5)
                ax1.axhline(a.tol, color="k", ls="--", lw=0.8); ax1.set_title("Resolved shear stress: max rel. error, full set", fontsize=10.5)
                k4 = ["T1", "T2", "T3", "T4"]; ax2.bar(k4 + ["push-fwd"], [ev[k] for k in k4] + [epf], color="C2"); ax2.set_yscale("log"); ax2.set_ylim(1e-18, 1)
                for i, v in enumerate([ev[k] for k in k4] + [epf]): ax2.text(i, v*2.2, f"{v:.1e}", ha="center", fontsize=8.5)
                ax2.axhline(a.tol, color="k", ls="--", lw=0.8); ax2.set_title("Current slip vectors $F_e s_\\alpha$, $F_e^{-T}n_\\alpha$ and\n$F_eL_pF_e^{-1}$ push-forward: max error", fontsize=10.5)
                fig.suptitle(f"Slip follows the lattice ({law} law). Lines: crystal starting at $\\alpha$ (slip systems rotated once). Markers: tests. "
                             "$V=\\mathrm{diag}(1.10,0.95)$, Cu, plane strain, $\\mu_c=754$ GPa", fontsize=11)
                fig.tight_layout(rect=(0, 0, 1, 0.96)); os.makedirs(a.fig, exist_ok=True)
                p = os.path.join(a.fig, f"slip_follows_lattice_{law}.png"); fig.savefig(p, dpi=130); print("  figure:", p)
            except ImportError:
                print("  (matplotlib not available: no figure)")
    print("\nRESULT:", "FAIL" if bad else "PASS"); sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
