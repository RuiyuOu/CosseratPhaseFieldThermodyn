#!/usr/bin/env python3
"""Slip-law comparison over a run: what resolved shear stress would MTS and PTW need for the slip rates the run produced,
compared with the stress the run actually carried (our slip law). Area-weighted max, 90th percentile and median over time.

For every snapshot k > 0 and every element:
  slip rate      g_a = |p_a(k) - p_a(k-1)| / (t_k - t_{k-1})       (p_a = accumulated signed slip, VTU cell data p1, p2, p3;
                                                                      a mean over the snapshot interval: reversals inside it cancel)
  dominant system a* = argmax_a g_a;  gdot = g_a*,  tau = |tau_a*(k)|  (VTU tau_1..tau_3, end of the interval)
  "active" element: gdot >= --rmin (default 1e-4 /ns = 1e5 /s)
Stress needed for that rate, per law (resolved shear stress, GPa):
  ours  inverse of 1/gdot = 1/(g0 x^m) + eta_eff/tau,  x = tau/(f_a tau_c),  eta_eff = eta_d sqrt(1 + (tau/(eta_d rho_m b v_c))^2)
        with the element's tau_c, eta_drag, rho_mob (VTU) and f_a = (sqrt3, 2/sqrt3, 2/sqrt3); infinite above the speed limit
        rho_m b v_c.  Compared with the actual |tau_a*| this checks the interval-mean rate against the law.
  MTS   Follansbee & Kocks (1988), OFHC Cu: sigma = sigma_a + S sigma_e mu(T)/mu0,
        S = [1 - (kT/(g0 mu b^3) ln(edot0/edot))^(1/q)]^(1/p), p = 2/3, q = 1; mu(T) = mu0 - D/(exp(T0/T) - 1).
        Valid up to edot0 = 1e7 /s; above, NaN (counted).  Threshold state (--mts-state):
          density (default)  sigma_e = M lambda_T mu0 b sqrt(rho_S + rho_G)   (Taylor, the element's own densities)
          fixed              sigma_e = --mts-sige (MPa)
  PTW   Preston, Tonks & Wallace (2003), Cu (--ptw preston) or the LANL half-hard Cu calibration (--ptw sjue), with the
        shear modulus G(T) = G0 (1 - alpha T/Tm) and xi_dot = (1/2)(4 pi n/3)^(1/3) sqrt(G/rho).  Strain state (--ptw-state):
          slip (default)  eps_p = slip_path_phys / M   (the element's own accumulated slip)
          fixed           eps_p = --ptw-ep
  MTS and PTW give a uniaxial flow stress sigma(edot); converted with a Taylor factor: tau = sigma/M, edot = gdot/M (--taylor 3.06).
  Temperature: the cell mean of the nodal T (VTU point data).
  Pressure (default on; --no-pressure for ambient): P = -(sig_11 + sig_22 + sig_33)/3 per element (Cauchy, GPa) and
  rho/rho0 = 1/J.  MTS: mu(P,T) = mu(T) + dG/dP P in the activation term and the mu/mu0 scaling.  PTW: G(P,T) =
  (G0 + dG/dP P)(1 - alpha T/Tm(rho)), Tm(rho) = Tm0 exp(2a(1 - rho0/rho)) (rho/rho0)^(2(g0 - a - 1/3)) (Lindemann,
  Steinberg-Guinan form), and xi_dot from the compressed rho and G.  dG/dP = 1.5 (--dGdP), g0 = 2.0, a = 1.5.
  No anchoring or rescaling: the published calibrations as they are.

Statistics per snapshot (area-weighted over active elements): max, p90, median of tau (actual), tau_ours, tau_MTS, tau_PTW,
gdot and gdot/(rho_m b v_c); ratios tau_MTS/tau and tau_PTW/tau; active area fraction; area fractions above the MTS validity
limit and above our speed limit.

usage:   python3 post/sliplaw_compare.py RUN_DIR [--out DIR] [--rmin 1e-4] [--taylor 3.06] [--mts-state density|fixed]
                 [--mts-sige 200] [--ptw preston|sjue] [--ptw-state slip|fixed] [--ptw-ep 0.1] [--no-pressure] [--dGdP 1.5]
                 [--every N] [--tmax T]
RUN_DIR holds <tag>.pvd and the <tag>_s*.vtu snapshots (and normally <tag>_meta.json, from which g0, m, B0, v_c, b, mu,
lambda_T are read when present).
Output (--out, default ./sliplaw_post): sliplaw_<tag>.csv (one row per snapshot), sliplaw_<tag>_summary.txt, sliplaw_<tag>.png.
Needs numpy, meshio, matplotlib (Python >= 3.8).  No scipy, no pandas.
"""
import argparse, glob, json, math, os, re, sys
import numpy as np
try:
    import meshio
except ImportError:
    sys.exit("sliplaw_compare.py needs meshio (pip install --user meshio)")

# ----------------------------------------------------------------------------- our law (defaults = production 7.5)
OURS = dict(g0=1e-6, m=50.0, B0=2e-5, vc=2.0, b=2.556e-4, mu=48.0, lamT=0.3)
FA = np.array([math.sqrt(3.0), 2.0/math.sqrt(3.0), 2.0/math.sqrt(3.0)])
# ----------------------------------------------------------------------------- MTS, OFHC Cu (Follansbee & Kocks 1988)
MTS = dict(sig_a=40.0, g0=1.6, edot0=1e7, p=2.0/3.0, q=1.0, mu0=51.3, D=3.0, T0=165.0, b=2.56e-10)
# ----------------------------------------------------------------------------- PTW, Cu
PTW_SETS = {
    "preston": dict(th=0.025, p=2.0, s0=0.0085, sinf=0.00055, kap=0.11, gam=1e-5, y0=1e-4, yinf=1e-4, y1=0.094, y2=0.575,
                    beta=0.25, G0=51.9, alpha=0.34, Tm=1356.0),
    "sjue":    dict(th=0.0283, p=1.0, s0=0.0047, sinf=0.002, kap=0.1191, gam=5e-7, y0=0.0024, yinf=0.0012, y1=0.0047,
                    y2=0.44, beta=0.44, G0=51.9, alpha=0.34, Tm=1356.0),
}
KB, NA = 1.380649e-23, 6.02214076e23
RHO_CU, MA_CU = 8960.0, 0.063546      # kg/m^3, kg/mol


def erf(x):
    """vectorised erf (Abramowitz-Stegun 7.1.26, |error| < 1.5e-7), no scipy"""
    x = np.asarray(x, float); s = np.sign(x); a = np.abs(x); t = 1.0/(1.0 + 0.3275911*a)
    y = 1.0 - (((((1.061405429*t - 1.453152027)*t) + 1.421413741)*t - 0.284496736)*t + 0.254829592)*t*np.exp(-a*a)
    return s*y


def ours_lograte(tau, tc, fa, etad, rhom, P):
    """ln |gdot| (1/ns) of the production law, overflow free"""
    x = np.maximum(tau, 1e-300)/(fa*tc)
    lth = math.log(P["g0"]) + P["m"]*np.log(x)
    etae = etad*np.sqrt(1.0 + (tau/(etad*rhom*P["b"]*P["vc"]))**2)
    ldr = np.log(np.maximum(tau, 1e-300)/etae)
    return -np.logaddexp(-lth, -ldr)


def ours_stress(gdot, tc, fa, etad, rhom, P, nit=90):
    """resolved shear stress (GPa) at which our law gives gdot; inf at or above the speed limit rho_m b v_c"""
    cap = rhom*P["b"]*P["vc"]; ok = gdot < cap*(1 - 1e-12)
    lo = np.full(gdot.shape, math.log(1e-8)); hi = np.full(gdot.shape, math.log(1e4))
    tgt = np.log(np.maximum(gdot, 1e-300))
    for _ in range(nit):
        mid = 0.5*(lo + hi)
        up = ours_lograte(np.exp(mid), tc, fa, etad, rhom, P) < tgt
        lo = np.where(up, mid, lo); hi = np.where(up, hi, mid)
    out = np.exp(0.5*(lo + hi)); out[~ok] = np.inf
    return out


def mts_sigma(edot, T, sig_e, Pr=0.0, dGdP=0.0):
    """MTS flow stress (MPa) at strain rate edot (1/s), temperature T (K), threshold sig_e (MPa) and pressure Pr (GPa);
    NaN above edot0.  mu(P,T) = mu(T) + dG/dP P enters both the activation term and the mu/mu0 scaling."""
    M = MTS
    mu = M["mu0"] - M["D"]/(np.exp(M["T0"]/T) - 1.0) + dGdP*Pr          # GPa
    mu = np.maximum(mu, 0.2*M["mu0"])
    X = KB*T/(M["g0"]*mu*1e9*M["b"]**3)*np.log(M["edot0"]/edot)
    S = np.clip(1.0 - np.clip(X, 0.0, None)**(1.0/M["q"]), 0.0, None)**(1.0/M["p"])
    s = M["sig_a"] + S*sig_e*mu/M["mu0"]
    return np.where(edot <= M["edot0"], s, np.nan)


def melt_T(Tm0, eta, g0=2.0, a=1.5):
    """Lindemann melt temperature vs compression eta = rho/rho0 (Steinberg-Guinan form):
    Tm = Tm0 exp(2a(1 - 1/eta)) eta^(2(g0 - a - 1/3))"""
    return Tm0*np.exp(2*a*(1.0 - 1.0/eta))*eta**(2*(g0 - a - 1.0/3.0))


def ptw_sigma(edot, T, ep, P, Pr=0.0, eta=1.0, dGdP=0.0, g0=2.0, a=1.5):
    """PTW flow stress (MPa) at strain rate edot (1/s), temperature T (K), plastic strain ep, pressure Pr (GPa) and
    compression eta = rho/rho0: G(P,T) = (G0 + dG/dP P)(1 - alpha T/Tm(eta)), Tm(eta) Lindemann, xi_dot from rho and G"""
    Tm = melt_T(P["Tm"], eta, g0, a)
    G = np.maximum((P["G0"] + dGdP*Pr)*(1.0 - P["alpha"]*T/Tm), 0.05*P["G0"])   # GPa
    rho = RHO_CU*eta; n = rho*NA/MA_CU
    xid = 0.5*(4*math.pi*n/3.0)**(1.0/3.0)*np.sqrt(G*1e9/rho)              # 1/s
    That = T/Tm; r = edot/(P["gam"]*xid); A = P["kap"]*That*np.log(P["gam"]*xid/edot)
    ts = np.maximum(P["s0"] - (P["s0"] - P["sinf"])*erf(A), P["s0"]*r**P["beta"])
    ty = np.maximum(P["y0"] - (P["y0"] - P["yinf"])*erf(A), np.minimum(P["y1"]*r**P["y2"], P["s0"]*r**P["beta"]))
    al = (P["s0"] - ty)/P["p"]
    with np.errstate(all="ignore"):
        be = (ts - ty)/al; ph = np.exp(be) - 1.0
        th = ts + al*np.log(1.0 - ph*np.exp(-be - P["th"]*ep/(al*ph)))
    th = np.where((ts <= ty + 1e-15) | ~np.isfinite(th), ts, th)
    return 2.0*th*G*1e3                                                   # MPa


def wq(v, w, q):
    """area-weighted quantile (q in [0,1]) of the finite entries"""
    m = np.isfinite(v)
    if not m.any():
        return float("nan")
    v, w = v[m], w[m]; i = np.argsort(v); v, w = v[i], w[i]; c = np.cumsum(w); c /= c[-1]
    return float(np.interp(q, c, v))


def read_pvd(run):
    pvds = sorted(glob.glob(os.path.join(run, "*.pvd")))
    if not pvds:
        sys.exit(f"no .pvd in {run}")
    pvd = pvds[0]; tag = os.path.basename(pvd)[:-4]
    ent = re.findall(r'timestep="([^"]+)"[^>]*file="([^"]+)"', open(pvd).read())
    best = {}
    for t, f in ent:
        t = round(float(t), 9); p = os.path.join(run, f)
        if not os.path.isfile(p):
            continue
        if t not in best or "_s" in os.path.basename(f):                  # prefer _s<k> over _t<k> at the same time
            best[t] = p
    return tag, sorted(best.items())


def load_params(run, tag):
    P = dict(OURS); mf = os.path.join(run, f"{tag}_meta.json")
    if os.path.isfile(mf):
        try:
            j = json.load(open(mf))
            P["g0"] = float(j["rate"]["gdot0_per_ns"]); P["m"] = float(j["rate"]["m"]); P["B0"] = float(j["rate"]["B0_GPa_ns"])
            h = j.get("hardening", {}); P["mu"] = float(h.get("mu", P["mu"])); P["b"] = float(h.get("b_um", P["b"]))
            P["lamT"] = float(h.get("lambda_T", P["lamT"]))
            md = j.get("mobile_drag", {})
            if "vc_um_per_ns" in md:
                P["vc"] = float(md["vc_um_per_ns"])
        except Exception as ex:
            print(f"  (meta.json not fully read: {ex}; defaults used)")
    return P


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("run"); ap.add_argument("--out", default="sliplaw_post")
    ap.add_argument("--rmin", type=float, default=1e-4, help="active if the dominant slip rate >= rmin (1/ns)")
    ap.add_argument("--taylor", type=float, default=3.06)
    ap.add_argument("--mts-state", default="density", choices=("density", "fixed")); ap.add_argument("--mts-sige", type=float, default=200.0)
    ap.add_argument("--ptw", default="preston", choices=tuple(PTW_SETS)); ap.add_argument("--ptw-state", default="slip", choices=("slip", "fixed"))
    ap.add_argument("--ptw-ep", type=float, default=0.1)
    ap.add_argument("--no-pressure", action="store_true", help="evaluate MTS and PTW at ambient pressure and density")
    ap.add_argument("--dGdP", type=float, default=1.5, help="pressure derivative of the shear modulus (Cu, Sjue 2023)")
    ap.add_argument("--lind-g0", type=float, default=2.0, help="Grueneisen gamma0 in the Lindemann melt law")
    ap.add_argument("--lind-a", type=float, default=1.5, help="first-order volume correction a in the Lindemann melt law")
    ap.add_argument("--every", type=int, default=1); ap.add_argument("--tmax", type=float, default=None)
    a = ap.parse_args()
    tag, snaps = read_pvd(a.run)
    if a.tmax is not None:
        snaps = [s for s in snaps if s[0] <= a.tmax + 1e-9]
    snaps = snaps[::a.every]
    if len(snaps) < 2:
        sys.exit("need at least two snapshots")
    P = load_params(a.run, tag); PT = PTW_SETS[a.ptw]; M = a.taylor
    os.makedirs(a.out, exist_ok=True)
    print(f"{tag}: {len(snaps)} snapshots, t = {snaps[0][0]} .. {snaps[-1][0]} ns; our law g0={P['g0']}/ns m={P['m']} "
          f"B0={P['B0']} v_c={P['vc']}; MTS state {a.mts_state}; PTW {a.ptw}, state {a.ptw_state}; Taylor M = {M}; "
          + ("ambient pressure" if a.no_pressure else f"pressure on (dG/dP = {a.dGdP}, Lindemann g0 = {a.lind_g0}, a = {a.lind_a})"))
    rows = []; prev = None; area = None; cells = None
    for k, (t, f) in enumerate(snaps):
        m = meshio.read(f)
        if area is None:
            cells = m.cells_dict["triangle"]; X = m.points[:, :2][cells]
            area = 0.5*np.abs((X[:, 1, 0] - X[:, 0, 0])*(X[:, 2, 1] - X[:, 0, 1]) - (X[:, 2, 0] - X[:, 0, 0])*(X[:, 1, 1] - X[:, 0, 1]))
        cd = {kk: np.asarray(v[0]).ravel() for kk, v in m.cell_data.items()}
        p = np.stack([cd["p1"], cd["p2"], cd["p3"]], 1)
        if prev is None:
            prev = (t, p); continue
        dt = t - prev[0]
        if dt <= 0:
            continue
        g = np.abs(p - prev[1])/dt; prev = (t, p)
        ia = np.argmax(g, 1); gd = g[np.arange(len(g)), ia]
        tau = np.abs(np.stack([cd["tau_1"], cd["tau_2"], cd["tau_3"]], 1)[np.arange(len(g)), ia])
        tc = cd["tau_c"]; etad = cd["eta_drag"]; rhom = cd["rho_mob"]; fa = FA[ia]
        Tn = np.asarray(m.point_data["T"]).ravel() if "T" in m.point_data else np.full(m.points.shape[0], 300.0)
        Tc = Tn[cells].mean(1)
        if a.no_pressure:
            Pr = np.zeros(len(Tc)); eta = np.ones(len(Tc))
        else:
            Pr = -(cd["sig_11"] + cd["sig_22"] + cd["sig_33"])/3.0           # GPa, compression positive
            eta = 1.0/np.clip(cd["J"], 0.3, 3.0)                            # rho/rho0 = 1/J
        act = gd >= a.rmin; w = area[act]; atot = area.sum()
        r = dict(t_ns=t, n_active=int(act.sum()), active_area_frac=float(w.sum()/atot))
        if act.any():
            g_a = gd[act]; ed = g_a*1e9/M                                           # macroscopic strain rate, 1/s
            t_act = tau[act]*1e3                                                    # MPa
            t_our = ours_stress(g_a, tc[act], fa[act], etad[act], rhom[act], P)*1e3
            rhoF = cd["rhoS"][act] + cd["rhoG"][act]
            sige = (M*P["lamT"]*MTS["mu0"]*1e3*P["b"]*np.sqrt(np.maximum(rhoF, 0.0))) if a.mts_state == "density" else np.full(ed.shape, a.mts_sige)
            dG = 0.0 if a.no_pressure else a.dGdP
            t_mts = mts_sigma(ed, Tc[act], sige, Pr[act], dG)/M
            ep = (cd["slip_path_phys"][act]/M) if a.ptw_state == "slip" else np.full(ed.shape, a.ptw_ep)
            t_ptw = ptw_sigma(ed, Tc[act], ep, PT, Pr[act], eta[act], dG, a.lind_g0, a.lind_a)/M
            cap = rhom[act]*P["b"]*P["vc"]
            for name, v in (("tau_actual", t_act), ("tau_ours", t_our), ("tau_MTS", t_mts), ("tau_PTW", t_ptw),
                            ("ratio_MTS", t_mts/t_act), ("ratio_PTW", t_ptw/t_act), ("ratio_ours", t_our/t_act),
                            ("gdot", g_a), ("speed_frac", g_a/cap)):
                r[f"{name}_max"] = float(np.nanmax(np.where(np.isfinite(v), v, np.nan))) if np.isfinite(v).any() else float("nan")
                r[f"{name}_p90"] = wq(v, w, 0.90); r[f"{name}_med"] = wq(v, w, 0.50)
            r["frac_above_MTS_limit"] = float(w[~np.isfinite(t_mts)].sum()/w.sum())
            r["frac_above_our_cap"] = float(w[~np.isfinite(t_our)].sum()/w.sum())
            r["T_med"] = wq(Tc[act], w, 0.5)
        rows.append(r)
        if k % 20 == 0 or k == len(snaps) - 1:
            print(f"  t = {t:7.3f} ns: active {r['active_area_frac']*100:5.1f} %"
                  + (f", median tau = {r['tau_actual_med']:.0f} MPa, MTS {r['tau_MTS_med']:.0f}, PTW {r['tau_PTW_med']:.0f}" if r["n_active"] else ""))
    keys = sorted({kk for r in rows for kk in r}, key=lambda s: (s != "t_ns", s))
    csvp = os.path.join(a.out, f"sliplaw_{tag}.csv")
    with open(csvp, "w") as fh:
        fh.write(",".join(keys) + "\n")
        for r in rows:
            fh.write(",".join(f"{r.get(kk, float('nan')):.6g}" for kk in keys) + "\n")
    # summary
    act_rows = [r for r in rows if r["n_active"]]
    with open(os.path.join(a.out, f"sliplaw_{tag}_summary.txt"), "w") as fh:
        fh.write(f"run {tag}; snapshots {len(rows)} (t = {rows[0]['t_ns']} .. {rows[-1]['t_ns']} ns); active: gdot >= {a.rmin}/ns\n")
        fh.write(f"MTS: Follansbee-Kocks OFHC Cu, state {a.mts_state}; PTW: {a.ptw}, state {a.ptw_state}; Taylor M = {M}; "
                 + ("ambient pressure\n" if a.no_pressure else f"pressure on: dG/dP = {a.dGdP}, Lindemann g0 = {a.lind_g0}, a = {a.lind_a}\n"))
        for q in ("tau_actual", "tau_ours", "tau_MTS", "tau_PTW", "ratio_MTS", "ratio_PTW", "gdot", "speed_frac"):
            v = np.array([r.get(q + "_med", np.nan) for r in act_rows]); vm = np.array([r.get(q + "_max", np.nan) for r in act_rows])
            if len(v):
                fh.write(f"{q:12s} median over time of the median {np.nanmedian(v):10.4g}; largest max {np.nanmax(vm):10.4g}\n")
        if act_rows:
            fh.write(f"largest area fraction above the MTS validity limit (1e7/s): {max(r['frac_above_MTS_limit'] for r in act_rows):.3f}\n")
            fh.write(f"largest area fraction at/above our speed limit: {max(r['frac_above_our_cap'] for r in act_rows):.3f}\n")
    # figure
    try:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available: no figure"); return
    from matplotlib.lines import Line2D
    T = np.array([r["t_ns"] for r in rows]); G = lambda key: np.array([r.get(key, np.nan) for r in rows])
    fig, ax = plt.subplots(2, 2, figsize=(15, 10), sharex=True); ax = ax.ravel()
    styles = [Line2D([], [], color="0.3", ls="-", lw=2, label="median"), Line2D([], [], color="0.3", ls="--", lw=1.3, label="90th percentile"),
              Line2D([], [], color="0.3", ls=":", lw=1.3, label="max")]
    def trio(x, key, c, lw=2.0):
        x.plot(T, G(key + "_med"), "-", color=c, lw=lw); x.plot(T, G(key + "_p90"), "--", color=c, lw=lw*0.65); x.plot(T, G(key + "_max"), ":", color=c, lw=lw*0.65)
    def legend(x, items, loc="upper right"):
        h = [Line2D([], [], color=c, lw=2.2, label=l) for c, l in items] + styles
        x.legend(handles=h, fontsize=8.5, loc=loc, ncol=2, framealpha=0.85)
    q0 = [("tau_actual", "k", "run (our law): actual $|\\tau|$"), ("tau_MTS", "C4", "MTS, Cu: needed $\\tau$"),
          ("tau_PTW", "C2", f"PTW, Cu ({a.ptw}): needed $\\tau$"), ("tau_ours", "0.6", "our law at the interval-mean rate")]
    for key, c, _ in q0:
        trio(ax[0], key, c, 2.2 if key != "tau_ours" else 1.3)
    legend(ax[0], [(c, l) for _, c, l in q0]); ax[0].set_yscale("log"); ax[0].set_ylabel("resolved shear stress (MPa)")
    ax[0].set_title("Stress carried by the run vs stress MTS / PTW need for the same slip rate", fontsize=11)
    q1 = [("ratio_MTS", "C4", "MTS / run"), ("ratio_PTW", "C2", "PTW / run"), ("ratio_ours", "0.6", "our law (interval rate) / run")]
    for key, c, _ in q1:
        trio(ax[1], key, c)
    ax[1].axhline(1, color="k", lw=0.8); ax[1].set_yscale("log"); ax[1].set_ylabel("needed / actual stress")
    legend(ax[1], [(c, l) for _, c, l in q1]); ax[1].set_title("Ratio of the needed stress to the actual stress (element by element)", fontsize=11)
    q2 = [("gdot", "C0", "slip rate $\\dot\\gamma$ (1/ns)"), ("speed_frac", "C3", "$\\dot\\gamma/(\\rho_m b v_c)$ (our speed limit = 1)")]
    for key, c, _ in q2:
        trio(ax[2], key, c)
    ax[2].axhline(1, color="C3", lw=0.7); ax[2].axhline(1e7*M*1e-9, color="C4", lw=0.9, ls="-.")
    ax[2].text(T.min(), 1e7*M*1e-9*1.15, " MTS validity limit ($10^7$/s macroscopic)", color="C4", fontsize=8.5)
    legend(ax[2], [(c, l) for _, c, l in q2], loc="lower right"); ax[2].set_yscale("log"); ax[2].set_ylabel("slip rate, speed fraction")
    ax[2].set_title("Slip rate of the dominant system (active elements)", fontsize=11)
    ax[3].plot(T, G("active_area_frac"), "k-", lw=2, label=f"active area fraction ($\\dot\\gamma\\geq${a.rmin:g}/ns)")
    ax[3].plot(T, G("frac_above_MTS_limit"), "-", color="C4", lw=1.6, label="of active: above the MTS validity limit")
    ax[3].plot(T, G("frac_above_our_cap"), "-", color="C3", lw=1.6, label="of active: at/above our speed limit")
    ax[3].set_ylim(-0.02, 1.02); ax[3].set_ylabel("area fraction"); ax[3].set_title("Coverage", fontsize=11); ax[3].legend(fontsize=8.5, loc="best")
    for x in ax:
        x.axvspan(0, 3, color="0.92", zorder=0); x.grid(alpha=0.3); x.set_xlim(0, max(T.max(), 0.1))
    for x in ax[2:]:
        x.set_xlabel("t (ns)  (shaded: loading, 0-3 ns)")
    fig.suptitle(f"{tag}: slip-law comparison (MTS state: {a.mts_state}, PTW: {a.ptw}/{a.ptw_state}, Taylor M = {M}; {'ambient pressure' if a.no_pressure else 'with pressure'}; unscaled)", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.96)); png = os.path.join(a.out, f"sliplaw_{tag}.png"); fig.savefig(png, dpi=130)
    print("wrote", csvp, "and", png)


if __name__ == "__main__":
    main()
