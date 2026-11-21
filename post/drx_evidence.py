#!/usr/bin/env python3
"""drx_evidence.py -- quantitative tests for dynamic recrystallization (DRX) in one run, against laws from experiment and DDD.

usage:  BASE=<tag> python3 drx_evidence.py <run_dir> <out_dir> [--every N] [--grains <deliver>/grains]
        (run_dir = out_<tag>; reads its VTU frames, <tag>_log.csv, <tag>_thermo.csv, <tag>_meta.json; incremental: frames already in
        <out_dir>/drx_grains.csv are skipped, FORCE=1 redoes all)
env:    GRAD (2 deg/um)  MERGE (2 deg)  MINAREA (0.1 um^2)  SYM (90)  HAGB (15 deg)  PHIIN (0.9)  NP (processes, 1)  FORCE=1

Segmentation: grain_count.py (plateaus of theta_lattice_deg, territories, merge), the same as grain_post.py. The mesh is Lagrangian, so the
cell field Label is the prior (parent) grain of every element at every time.

Per grain and frame (drx_grains.csv): territory area and equivalent diameter, parent, kind (remnant = the largest grain descended from its
parent, new = every other descendant), orientation spread GOS over the interior (territory elements with phi >= PHIIN: like an EBSD grain,
it includes smooth orientation ramps) and over the plateau only (gos_plateau_deg), and interior means
of rho_S, rho_G, slip_path_phys, wstar, |theta - theta0_parent|, T, tau_c, von Mises stress.  Per boundary (drx_boundaries.csv):
misorientation (raw, mod SYM), length, intra (both grains from the same parent: deformation-induced) or prior.

Figures (out_dir):
  drx1_defects.png       new vs remnant grains: area-weighted rho_S + rho_G and GOS against time   (recrystallized = clean)
  drx2_derby.png         D/b against sigma/mu of new grains, Derby (1991) band K = 1..10           (sigma = sqrt(3) tau_c)
  drx3_similitude.png    grain size against (rho_S + rho_G)^-1/2, lines d = 5, 10, 20 rho^-1/2     (similitude, DDD/experiment)
  drx4_misorientation.png deformation-induced boundaries: theta_av p(theta/theta_av) per frame (Hughes et al. 1997 collapse) and
                         theta_av against the mean slip path (slope 1/2 reference)
  drx5_hagb.png          length fraction of deformation-induced boundaries >= HAGB and number of new grains against time
  drx6_sites.png         distance of new-grain centres to the prior boundaries, against uniformly placed points  (nucleation sites)
  drx7_energy.png        stored energy, dissipation and stored fraction E_stored / (E_stored + D_cum) against time
  drx8_mechanism.png     new grains: |w*| (grain-boundary spin) against |theta - theta0_parent| (total lattice rotation)
  drx_summary.json       the numbers behind the figures
"""
import sys, os, re, json, glob, argparse, csv
from multiprocessing import Pool
import numpy as np
from scipy.spatial import cKDTree
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import grain_count as gc
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

E = os.environ.get
PAR = argparse.Namespace(grad=float(E("GRAD", 2)), merge=float(E("MERGE", 2)), minarea=float(E("MINAREA", 0.1)), sym=float(E("SYM", 90)),
                         lo=2.0, mid=5.0, hi=15.0, scan=False)
HAGB, PHIIN = float(E("HAGB", 15)), float(E("PHIIN", 0.9))
GCOLS = ["t_ns", "step", "id", "parent", "kind", "area_um2", "d_um", "cx_um", "cy_um", "theta_deg", "rot_from_parent_deg", "gos_deg", "gos_plateau_deg",
         "rhoS", "rhoG", "slip", "wstar_deg", "T", "tau_c", "sig_eq"]
BCOLS = ["t_ns", "step", "i", "j", "intra", "misorient_deg", "misorient_sym_deg", "length_um"]


def frames(run, base):
    """(t, file) of every VTU listed in <base>.pvd, plus the _t snapshots, sorted by time."""
    out = {}
    pvd = os.path.join(run, base + ".pvd")
    if os.path.exists(pvd):
        for t, f in re.findall(r'timestep="([^"]+)"[^>]*file="([^"]+)"', open(pvd).read()):
            p = os.path.join(run, os.path.basename(f))
            if os.path.exists(p): out[round(float(t), 6)] = p
    for p in glob.glob(os.path.join(run, base + "_[st]*.vtu")):
        m = re.search(r"_([st])(\d+)\.vtu$", p)
        if m: out.setdefault(round(int(m.group(2)) * 0.01, 6), p)
    return sorted(out.items())


def analyse(args):
    t, fn, theta0 = args
    d = gc.read_vtu(fn)
    pts = d["points"][:, :2].astype(float)
    conn, off, typ = d["connectivity"], d["offsets"], d["types"]
    start = np.r_[0, off[:-1]]; istri = typ == 5
    tri = np.array([conn[s:s + 3] for s, ty in zip(start, typ) if ty == 5])
    th, phi = d["theta_lattice_deg"].astype(float).ravel(), d["phi"].astype(float).ravel()
    M = gc.Mesh(pts, tri, th, phi)
    seg = gc.segment(M, PAR.grad, PAR.minarea, PAR.merge)
    if seg is None: return t, fn, [], []
    c = lambda k: np.asarray(d[k], float).ravel()[istri]
    lab = d["Label"][istri].astype(int)
    rho_s, rho_g, slip, wst, tauc = c("rhoS"), c("rhoG"), c("slip_path_phys"), c("wstar"), c("tau_c")
    s11, s22, s33, s12 = c("sig_11"), c("sig_22"), c("sig_33"), 0.5 * (c("sig_12") + c("sig_21"))
    seq = np.sqrt(0.5 * ((s11 - s22) ** 2 + (s22 - s33) ** 2 + (s33 - s11) ** 2) + 3 * s12 ** 2)
    Te, phie = d["T"].astype(float).ravel()[tri].mean(1), phi[tri].mean(1)
    terr, mu, A = seg["terr"], seg["mu"], M.A
    G = seg["G"]
    par = np.array([np.bincount(lab[terr == g], weights=A[terr == g], minlength=lab.max() + 1).argmax() if (terr == g).any() else -1 for g in range(G)])
    area = np.array([A[terr == g].sum() for g in range(G)])
    kind = np.array(["new"] * G, dtype=object)
    for p in set(par):
        ids = np.flatnonzero(par == p); kind[ids[np.argmax(area[ids])]] = "remnant"
    step = int(re.search(r"_[st](\d+)\.vtu$", fn).group(1)) if re.search(r"_[st](\d+)\.vtu$", fn) else -1
    rows = []
    for g in range(G):
        t_ = terr == g; inn = t_ & (phie >= PHIIN)
        if inn.sum() < 3: inn = t_
        w = A[inn]; m = lambda x: float((x[inn] * w).sum() / w.sum())
        gos = float(np.sqrt((w * (M.tm[inn] - m(M.tm)) ** 2).sum() / w.sum()))
        pl = seg["plat"] == g; wp = A[pl]
        gosp = float(np.sqrt((wp * (M.tm[pl] - (M.tm[pl] * wp).sum() / wp.sum()) ** 2).sum() / wp.sum())) if pl.any() else np.nan
        cen = (M.cen[t_] * A[t_, None]).sum(0) / A[t_].sum()
        th0 = np.degrees(theta0[par[g] - 1]) if 1 <= par[g] <= len(theta0) else np.nan
        rows.append(dict(t_ns=t, step=step, id=g, parent=int(par[g]), kind=kind[g], area_um2=area[g], d_um=2 * np.sqrt(area[g] / np.pi),
                         cx_um=cen[0], cy_um=cen[1], theta_deg=mu[g], rot_from_parent_deg=abs(mu[g] - th0), gos_deg=gos, gos_plateau_deg=gosp,
                         rhoS=m(rho_s), rhoG=m(rho_g), slip=m(slip), wstar_deg=np.degrees(abs(m(wst))), T=m(Te), tau_c=m(tauc), sig_eq=m(seq)))
    _, bnd = gc.analyse(M, seg, lab, PAR.sym)
    brows = [dict(t_ns=t, step=step, i=b["i"], j=b["j"], intra=int(par[b["i"]] == par[b["j"]]), misorient_deg=b["misorient_deg"],
                  misorient_sym_deg=b["misorient_sym_deg"], length_um=b["length_um"]) for b in bnd]
    return t, fn, rows, brows


def prior_boundary_points(fn):
    """Midpoints of the mesh edges between elements of different Label (the prior grain boundaries), reference coordinates."""
    d = gc.read_vtu(fn); pts = d["points"][:, :2].astype(float)
    conn, off, typ = d["connectivity"], d["offsets"], d["types"]; start = np.r_[0, off[:-1]]
    tri = np.array([conn[s:s + 3] for s, ty in zip(start, typ) if ty == 5]); lab = d["Label"][typ == 5].astype(int)
    M = gc.Mesh(pts, tri, np.zeros(len(pts)), np.zeros(len(pts)))
    dif = lab[M.pa] != lab[M.pb]
    return M.emid[dif], pts.min(0), pts.max(0)


def read_csv(fn):
    if not os.path.exists(fn): return []
    with open(fn) as f: return list(csv.DictReader(f))


def num(rows, k): return np.array([float(r[k]) for r in rows])


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("run"); ap.add_argument("out"); ap.add_argument("--every", type=int, default=1)
    ap.add_argument("--grains", default=None, help="deliver_<tag>/grains (for the grain count at every frame)")
    a = ap.parse_args()
    base = E("BASE") or re.sub(r"^out_", "", os.path.basename(os.path.normpath(a.run)))
    os.makedirs(a.out, exist_ok=True)
    meta = json.load(open(os.path.join(a.run, base + "_meta.json")))
    theta0 = np.array(meta["grain_theta0_rad"], float); mu_s, b = meta["hardening"]["mu"], meta["hardening"]["b_um"]
    gf, bf = os.path.join(a.out, "drx_grains.csv"), os.path.join(a.out, "drx_boundaries.csv")
    if E("FORCE") == "1":
        for f in (gf, bf):
            if os.path.exists(f): os.remove(f)
    done = {round(float(r["t_ns"]), 6) for r in read_csv(gf)}
    fr = frames(a.run, base)[::a.every]
    todo = [(t, f, theta0) for t, f in fr if round(t, 6) not in done]
    print(f"{len(fr)} frames, {len(todo)} to analyse")
    res = Pool(int(E("NP", 1))).map(analyse, todo) if int(E("NP", 1)) > 1 else map(analyse, todo)
    for t, fn, rows, brows in res:
        for f, cols, rr in ((gf, GCOLS, rows), (bf, BCOLS, brows)):
            new = not os.path.exists(f)
            with open(f, "a", newline="") as fh:
                w = csv.DictWriter(fh, cols); new and w.writeheader(); w.writerows(rr)
        print(f"  t = {t:6.2f} ns: {len(rows)} grains ({sum(r['kind'] == 'new' for r in rows)} new)", flush=True)
    G, B = read_csv(gf), read_csv(bf)
    if not G: print("no grains"); return
    S = {"run": base, "params": dict(vars(PAR), hagb=HAGB, phi_interior=PHIIN), "mu_GPa": mu_s, "b_um": b}
    t = num(G, "t_ns"); ts = np.unique(t); new = np.array([r["kind"] == "new" for r in G])
    area, d_um, rho, gos = num(G, "area_um2"), num(G, "d_um"), num(G, "rhoS") + num(G, "rhoG"), num(G, "gos_deg")
    wavg = lambda x, m: float((x[m] * area[m]).sum() / area[m].sum()) if m.any() else np.nan

    # 1. defect content and orientation spread
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    for lab_, msk, col in (("new grains", new, "C3"), ("remnants of the prior grains", ~new, "C0")):
        ax[0].plot(ts, [wavg(rho, msk & (t == x)) for x in ts], "o-", color=col, label=lab_)
        ax[1].plot(ts, [wavg(gos, msk & (t == x)) for x in ts], "o-", color=col, label=lab_)
    ax[0].set_ylabel("area-weighted ρ_S + ρ_G (µm⁻²)"); ax[1].set_ylabel("area-weighted GOS (deg)")
    ax[1].axhline(2, color="k", ls=":", lw=1); ax[1].text(ts[0], 2.05, "EBSD threshold ≈ 1–2°", fontsize=8)
    for x in ax: x.set_xlabel("t (ns)"); x.legend(fontsize=8); x.grid(alpha=.3)
    fig.suptitle(f"{base}: are the new grains clean?"); fig.tight_layout(); fig.savefig(os.path.join(a.out, "drx1_defects.png"), dpi=120); plt.close(fig)
    last = t == ts[-1]
    S["defects_last"] = dict(new=dict(n=int((new & last).sum()), rho=wavg(rho, new & last), gos=wavg(gos, new & last)),
                             remnant=dict(n=int((~new & last).sum()), rho=wavg(rho, ~new & last), gos=wavg(gos, ~new & last)))

    # 2. Derby: D/b vs sigma/mu
    sig = np.sqrt(3) * num(G, "tau_c")
    fig, ax = plt.subplots(figsize=(5.5, 4.5))
    xs = np.logspace(-4, -1, 50)
    for K in (1, 10): ax.plot((K / xs) ** 1.5, xs, "k--", lw=1)
    sc = ax.scatter(d_um[new] / b, sig[new] / mu_s, c=t[new], s=14, cmap="viridis", label="new grains")
    ax.scatter(d_um[~new & last] / b, sig[~new & last] / mu_s, marker="s", facecolor="none", edgecolor="C0", s=22, label="remnants (last frame)")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("D / b"); ax.set_ylabel("σ / μ   (σ = √3 τ_c)")
    ax.text(0.03, 0.04, "dashed: Derby (1991), (σ/μ)(D/b)^{2/3} = K, K = 1 and 10", transform=ax.transAxes, fontsize=7)
    fig.colorbar(sc, label="t (ns)"); ax.legend(fontsize=7, loc="upper right"); ax.grid(alpha=.3, which="both")
    fig.tight_layout(); fig.savefig(os.path.join(a.out, "drx2_derby.png"), dpi=120); plt.close(fig)
    Kd = (sig / mu_s) * (d_um / b) ** (2 / 3)
    S["derby"] = dict(K_median_new=float(np.median(Kd[new])) if new.any() else None,
                      K_quartiles_new=[float(np.percentile(Kd[new], q)) for q in (25, 75)] if new.any() else None,
                      frac_new_in_band=float(np.mean((Kd[new] >= 1) & (Kd[new] <= 10))) if new.any() else None)

    # 3. similitude: d vs rho^-1/2
    fig, ax = plt.subplots(figsize=(5.5, 4.5))
    ls = rho ** -0.5
    ax.scatter(ls[new], d_um[new], c=t[new], s=14, cmap="viridis", label="new grains")
    ax.scatter(ls[~new & last], d_um[~new & last], marker="s", facecolor="none", edgecolor="C0", s=22, label="remnants (last frame)")
    xx = np.logspace(np.log10(max(ls.min(), 1e-3) / 2), np.log10(ls.max() * 2), 20)
    for D_, st in ((5, ":"), (10, "--"), (20, ":")): ax.plot(xx, D_ * xx, "k" + st, lw=1); ax.text(xx[-1], D_ * xx[-1], f"{D_}", fontsize=7)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("(ρ_S + ρ_G)^{-1/2} (µm)"); ax.set_ylabel("grain diameter d (µm)")
    ax.set_title("similitude: d = D ρ^{-1/2}, D ≈ 10", fontsize=9); ax.legend(fontsize=7); ax.grid(alpha=.3, which="both")
    fig.tight_layout(); fig.savefig(os.path.join(a.out, "drx3_similitude.png"), dpi=120); plt.close(fig)
    S["similitude"] = dict(D_median_new=float(np.median(d_um[new] / ls[new])) if new.any() else None)

    # 4. misorientation of deformation-induced boundaries: Hughes scaling and theta_av vs slip
    tb, intra, mis, L = num(B, "t_ns"), num(B, "intra") > 0, num(B, "misorient_sym_deg"), num(B, "length_um")
    slip_t = np.array([wavg(num(G, "slip"), t == x) for x in ts])
    fig, ax = plt.subplots(1, 2, figsize=(11, 4)); thav = []
    for k, x in enumerate(ts):
        m = intra & (tb == x) & (mis >= 0.5)
        if m.sum() < 3: thav.append(np.nan); continue
        av = float((mis[m] * L[m]).sum() / L[m].sum()); thav.append(av)
        h, e = np.histogram(mis[m] / av, bins=np.linspace(0, 4, 13), weights=L[m], density=True)
        ax[0].plot(0.5 * (e[1:] + e[:-1]), h, "-", color=plt.cm.viridis(k / max(len(ts) - 1, 1)), lw=1.2, label=f"{x:g} ns" if k % max(1, len(ts) // 6) == 0 else None)
    ax[0].set_xlabel("θ / θ_av"); ax[0].set_ylabel("θ_av p(θ)  (length-weighted)"); ax[0].legend(fontsize=7); ax[0].grid(alpha=.3)
    ax[0].set_title("deformation-induced boundaries (same parent): collapse?", fontsize=9)
    thav = np.array(thav); ok = np.isfinite(thav) & (slip_t > 0)
    ax[1].loglog(slip_t[ok], thav[ok], "o-", label="θ_av")
    if ok.sum() >= 2:
        x0 = slip_t[ok][0]; ax[1].loglog(slip_t[ok], thav[ok][0] * (slip_t[ok] / x0) ** 0.5, "k--", lw=1, label="slope 1/2")
        S["theta_av_slope"] = float(np.polyfit(np.log(slip_t[ok]), np.log(thav[ok]), 1)[0])
    ax[1].set_xlabel("area-mean slip path"); ax[1].set_ylabel("θ_av (deg)"); ax[1].legend(fontsize=8); ax[1].grid(alpha=.3, which="both")
    fig.tight_layout(); fig.savefig(os.path.join(a.out, "drx4_misorientation.png"), dpi=120); plt.close(fig)

    # 5. HAGB fraction of deformation-induced boundaries, number of new grains
    fh = [float(L[intra & (tb == x) & (mis >= HAGB)].sum() / L[intra & (tb == x)].sum()) if (intra & (tb == x)).any() else 0 for x in ts]
    nn = [int((new & (t == x)).sum()) for x in ts]
    fig, ax = plt.subplots(figsize=(6.5, 4)); ax.plot(ts, fh, "o-", color="C3"); ax.set_ylabel(f"length fraction ≥ {HAGB:g}° (deformation-induced)", color="C3")
    ax2 = ax.twinx(); ax2.plot(ts, nn, "s-", color="C0"); ax2.set_ylabel("new grains", color="C0")
    if a.grains and os.path.exists(os.path.join(a.grains, "grain_history.csv")):
        h = read_csv(os.path.join(a.grains, "grain_history.csv")); ax2.plot(num(h, "t_ns"), num(h, "n_grains") - len(theta0), "-", color="C0", alpha=.4, lw=1, label="grains − 8 (every frame)"); ax2.legend(fontsize=7)
    ax.set_xlabel("t (ns)"); ax.grid(alpha=.3); fig.tight_layout(); fig.savefig(os.path.join(a.out, "drx5_hagb.png"), dpi=120); plt.close(fig)
    S["hagb_fraction_intra"] = dict(zip([f"{x:g}" for x in ts], fh)); S["new_grains"] = dict(zip([f"{x:g}" for x in ts], nn))

    # 6. nucleation sites: distance to the prior boundaries vs uniform points
    pb, lo, hi = prior_boundary_points(fr[0][1]); tree = cKDTree(pb)
    cxy = np.c_[num(G, "cx_um"), num(G, "cy_um")]
    dn = tree.query(cxy[new])[0] if new.any() else np.array([])
    rng = np.random.default_rng(0); ru = tree.query(lo + (hi - lo) * rng.random((20000, 2)))[0]
    fig, ax = plt.subplots(figsize=(5.5, 4))
    ax.plot(np.sort(ru), np.linspace(0, 1, len(ru)), "k--", label="uniformly placed points")
    if len(dn): ax.plot(np.sort(dn), np.linspace(0, 1, len(dn)), "C3-", label=f"new-grain centres (n = {len(dn)})")
    ax.set_xlabel("distance to the nearest prior grain boundary (µm)"); ax.set_ylabel("cumulative fraction"); ax.legend(fontsize=8); ax.grid(alpha=.3)
    fig.tight_layout(); fig.savefig(os.path.join(a.out, "drx6_sites.png"), dpi=120); plt.close(fig)
    S["sites"] = dict(median_new_um=float(np.median(dn)) if len(dn) else None, median_uniform_um=float(np.median(ru)))
    if len(dn) >= 5:
        from scipy.stats import ks_2samp
        S["sites"]["ks_p"] = float(ks_2samp(dn, ru).pvalue)

    # 7. energetics
    lg = read_csv(os.path.join(a.run, base + "_log.csv"))
    if lg:
        tt, Es, Dc, Wp = num(lg, "t_ns"), num(lg, "E_stored"), num(lg, "D_cum"), num(lg, "W_plate")
        fig, ax = plt.subplots(1, 2, figsize=(11, 4))
        ax[0].plot(tt, Es, label="E_stored"); ax[0].plot(tt, Dc, label="D_cum (dissipated)"); ax[0].plot(tt, Wp, label="W_plate")
        ax[0].set_yscale("symlog", linthresh=1e-3); ax[0].set_xlabel("t (ns)"); ax[0].set_ylabel("GPa µm²"); ax[0].legend(fontsize=8); ax[0].grid(alpha=.3)
        dE = Es - Es[0]; fr_ = np.where(dE + Dc > 0, dE / np.maximum(dE + Dc, 1e-30), np.nan)
        ax[1].plot(tt, fr_, "C3"); ax[1].set_xlabel("t (ns)"); ax[1].set_ylabel("(E_stored − E_stored,0) / (that + D_cum)"); ax[1].grid(alpha=.3)
        ax[1].set_title("stored fraction of the plastic work (1 − β)", fontsize=9)
        fig.tight_layout(); fig.savefig(os.path.join(a.out, "drx7_energy.png"), dpi=120); plt.close(fig)
        k = int(np.nanargmax(Es)); S["energy"] = dict(E_stored_peak=float(Es[k]), t_peak_ns=float(tt[k]), E_stored_end=float(Es[-1]),
                                                       drop_after_peak=float(1 - Es[-1] / Es[k]) if Es[k] > 0 else None)

    # 8. mechanism: grain-boundary spin vs total lattice rotation of the new grains
    ws, rot = num(G, "wstar_deg"), num(G, "rot_from_parent_deg")
    fig, ax = plt.subplots(figsize=(5, 4.5))
    ax.scatter(rot[new], ws[new], c=t[new], s=14, cmap="viridis"); mx = max(rot[new].max() if new.any() else 1, 1)
    ax.plot([0, mx], [0, mx], "k--", lw=1); ax.set_xlabel("|θ − θ0,parent| (deg)"); ax.set_ylabel("|w*| (deg)")
    ax.set_title("new grains: rotation carried by the grain-boundary spin", fontsize=9); ax.grid(alpha=.3)
    fig.tight_layout(); fig.savefig(os.path.join(a.out, "drx8_mechanism.png"), dpi=120); plt.close(fig)
    if new.any(): S["mechanism"] = dict(median_wstar_over_rotation=float(np.median(ws[new] / np.maximum(rot[new], 1e-6))))
    json.dump(S, open(os.path.join(a.out, "drx_summary.json"), "w"), indent=1, default=float)
    print(json.dumps(S, indent=1, default=float))


if __name__ == "__main__":
    main()
