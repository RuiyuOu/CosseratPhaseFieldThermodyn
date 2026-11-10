#!/usr/bin/env python3
"""Resolved shear stress statistics from the saved VTU series: tau_1, tau_2, tau_3 and tau_c.

tau_a (a = 1,2,3) is the Mandel resolved shear stress M:(s_a x n_a) on slip system a, one value per triangle (GPa, signed);
tau_c is the slip strength of the element (GPa, includes hardening and thermal softening).  Both are P0 cell data written by
the solver in every <out>_<tag>.vtu.  Everything below is area weighted (the triangles have different areas).

What is computed, for every snapshot and every selection of elements:
  tau_max  = max_a |tau_a|  (largest resolved shear stress of the element)      -> max, p99, p95, p90, median, mean
  per system a: max |tau_a|, p99, mean |tau_a|, and the signed minimum / maximum
  tau_c                                                                         -> min, median, mean, max
  ratio tau_max / tau_c                                                         -> median, p99, area fraction >= 1
  area fraction with tau_max above the thresholds --thr (default 0.1 0.2 0.5 GPa; 0.5 GPa is the DDD value)
  temperature (cell mean of T): max, p99, median, mean;  pressure -(s11+s22+s33)/3: p01, median, p99, max, min
  summary: peaks, and the time/rate of the largest rise of mean T and of max T between snapshots (heating jump)
  mismatch beta = phi_e - theta (Cosserat-penalty residual) from u, Fp and dtheta_deg: max, p99, median (deg), and max|dtheta|;
    --log CSV overlays/summarises max_abs_dtheta_deg from the solver log
  coverage check: the expected end time (--tend, default 30 ns) and the largest gap between snapshots are reported
Selections:  all elements;  bulk (cell-mean phi >= --phisplit);  core (cell-mean phi < --phisplit);  optional --xrange.

Time of a snapshot: from the <out>.pvd next to the file when it lists it, otherwise from the file name
(_s<n> = loaded step n -> n*dt;  _t<k> = k/100 ns;  _t0 = 0;  other tags are skipped unless --all-tags).

usage:
  python3 tau_stats.py RUN_DIR [RUN_DIR2 ...] [--out DIR] [--dt 0.01] [--window 0.5 2.5] [--phisplit 0.5]
                       [--thr 0.1 0.2 0.5] [--xrange X0 X1] [--every N] [--tmax T] [--all-tags] [--no-plots]
  python3 tau_stats.py path/to/file.vtu            (a single snapshot)
RUN_DIR holds the <tag>_s*.vtu files (and normally <tag>.pvd).  Several run directories are compared in one figure.
Output (in --out, default ./tau_post): tau_stats_<tag>.csv (one row per snapshot and selection), tau_summary.csv (one row per run:
peaks and window means), tau_vs_time.png and, for several runs, tau_compare.png.
Needs numpy, meshio and (for the figures) matplotlib; no FreeFem, no pandas.
"""
import argparse, csv, glob, os, re, sys
import numpy as np

try:
    import meshio
except ImportError:
    sys.exit("tau_stats.py needs meshio (pip install --user meshio)")


# ------------------------------------------------------------------ small helpers
def wquantile(x, w, qs):
    """area-weighted quantiles (qs in [0,1]) of x"""
    if x.size == 0:
        return [np.nan] * len(qs)
    o = np.argsort(x)
    xs, ws = x[o], w[o]
    c = np.cumsum(ws) - 0.5 * ws
    c /= ws.sum()
    return list(np.interp(qs, c, xs))


def wmean(x, w):
    return float(np.sum(x * w) / np.sum(w)) if x.size else np.nan


def snapshot_time(path, pvd_times, dt, all_tags):
    """(time in ns or None, tag)"""
    base = os.path.basename(path)
    if base in pvd_times:
        m = re.search(r"_([st]\d+)\.vtu$", base)
        return pvd_times[base], (m.group(1) if m else base)
    m = re.search(r"_s(\d+)\.vtu$", base)
    if m:
        return int(m.group(1)) * dt, "s" + m.group(1)
    m = re.search(r"_t(\d+)\.vtu$", base)
    if m:
        return int(m.group(1)) / 100.0, "t" + m.group(1)
    m = re.search(r"_([A-Za-z]+\d*)\.vtu$", base)
    if m and all_tags:
        return None, m.group(1)
    return None, None


def read_pvd(directory):
    """file name -> time (ns) from every .pvd in the directory"""
    out = {}
    for pvd in glob.glob(os.path.join(directory, "*.pvd")):
        try:
            txt = open(pvd).read()
        except IOError:
            continue
        for tm, fn in re.findall(r'timestep="([^"]+)"\s+file="([^"]+)"', txt):
            try:
                out[os.path.basename(fn)] = float(tm)
            except ValueError:
                pass
    return out


def heating_jump(t, y):
    """(time, rate) of the largest increase rate between consecutive snapshots; time = middle of that interval"""
    ok = np.isfinite(t) & np.isfinite(y)
    t, y = t[ok], y[ok]
    if t.size < 2:
        return None
    dt = np.diff(t)
    good = dt > 0
    if not good.any():
        return None
    rate = np.where(good, np.diff(y) / np.where(good, dt, 1.0), -np.inf)
    i = int(np.argmax(rate))
    return float(0.5 * (t[i] + t[i + 1])), float(rate[i])


def coverage_report(tag, rows, tend):
    t = np.array(sorted(set(r["t_ns"] for r in rows if not np.isnan(r["t_ns"]))))
    if t.size == 0:
        return
    gaps = np.diff(t) if t.size > 1 else np.array([0.0])
    msg = "  coverage %s: %d snapshots, t = %.3f ... %.3f ns, largest gap %.3f ns" % (tag, t.size, t[0], t[-1], gaps.max())
    sys.stderr.write(msg + "\n")
    if tend is not None and t[-1] < tend - max(1.5 * gaps.max(), 0.05):
        sys.stderr.write("  WARNING: the series ends at %.3f ns, before the expected %.3f ns (missing files or a stopped run)\n" % (t[-1], tend))
    if gaps.max() > 0.5:
        sys.stderr.write("  WARNING: gap of %.2f ns between snapshots; a short heating jump could fall between them (-vtkevery smaller)\n" % gaps.max())


def read_log(path):
    """solver log CSV: returns dict with t_ns and max_abs_dtheta_deg (and a few others) as arrays, or None"""
    try:
        with open(path) as fh:
            head = fh.readline().strip().split(",")
        if "t_ns" not in head or "max_abs_dtheta_deg" not in head:
            sys.stderr.write("  log %s has no t_ns / max_abs_dtheta_deg column\n" % path)
            return None
        d = np.genfromtxt(path, delimiter=",", names=True, dtype=float, invalid_raise=False)
        return {"t_ns": np.atleast_1d(d["t_ns"]), "max_abs_dtheta_deg": np.atleast_1d(d["max_abs_dtheta_deg"])}
    except Exception as e:
        sys.stderr.write("  cannot read log %s: %s\n" % (path, e))
        return None


def run_tag(directory, files):
    base = os.path.basename(files[0])
    return re.sub(r"_(s\d+|t\d+|relaxed|stall\d+|tstop\d+)\.vtu$", "", base)


def run_label(tag):
    """velocity (m/s) and plate angle parsed from the tag, e.g. ..._a90_up400ms"""
    v = re.search(r"up(\d+)ms", tag)
    a = re.search(r"_a(\d+)(?:_|$)", tag)
    return (int(v.group(1)) if v else None, int(a.group(1)) if a else None)


# ------------------------------------------------------------------ one snapshot
COLS_SEL = ["all", "bulk", "core"]


def load_snapshot(path):
    m = meshio.read(path)
    if "triangle" not in m.cells_dict:
        raise ValueError("no triangle cells")
    tri = m.cells_dict["triangle"]
    cd = m.cell_data_dict
    for k in ("tau_1", "tau_2", "tau_3", "tau_c"):
        if k not in cd or "triangle" not in cd[k]:
            raise ValueError("missing cell field " + k)
    p = m.points[:, :2]
    x, y = p[tri, 0], p[tri, 1]
    area = 0.5 * np.abs((x[:, 1] - x[:, 0]) * (y[:, 2] - y[:, 0]) - (x[:, 2] - x[:, 0]) * (y[:, 1] - y[:, 0]))
    nt = len(tri)
    col = lambda name: np.asarray(cd[name]["triangle"], float).reshape(nt, -1)[:, 0]     # (nt,) or (nt,1) -> (nt,)
    tau = np.vstack([col("tau_%d" % i) for i in (1, 2, 3)])
    tauc = col("tau_c")
    if "phi" in m.point_data:
        phip = np.asarray(m.point_data["phi"], float).reshape(len(p), -1)[:, 0]     # (N,) or (N,1) -> (N,)
        phi = phip[tri].mean(1)
    else:
        phi = np.ones(len(tri))
    xc = x.mean(1)
    if "T" in m.point_data:
        Tp = np.asarray(m.point_data["T"], float).reshape(len(p), -1)[:, 0]
        Tc = Tp[tri].mean(1)
    else:
        Tc = None
    if all(k in cd and "triangle" in cd[k] for k in ("sig_11", "sig_22", "sig_33")):
        press = -(col("sig_11") + col("sig_22") + col("sig_33")) / 3.0       # GPa, compression positive
    else:
        press = None
    beta = dth = None
    need = ("Fp_11", "Fp_12", "Fp_21", "Fp_22")
    if "u" in m.point_data and "dtheta_deg" in m.point_data and all(k in cd and "triangle" in cd[k] for k in need):
        # F = I + grad u (P1 gradient per triangle), G = F Fp^-1, phi_e = atan2(G21-G12, G11+G22)
        # beta = phi_e - dtheta, wrapped to (-180,180]; this is the Cosserat-penalty mismatch
        u = np.asarray(m.point_data["u"], float).reshape(len(p), -1)[:, :2]
        dU = np.stack([u[tri][:, 1] - u[tri][:, 0], u[tri][:, 2] - u[tri][:, 0]], axis=2)
        F = np.eye(2) + dU @ np.linalg.inv(np.stack([np.stack([x[:, 1] - x[:, 0], y[:, 1] - y[:, 0]], 1),
                                                     np.stack([x[:, 2] - x[:, 0], y[:, 2] - y[:, 0]], 1)], axis=2))
        Fp = np.stack([np.stack([col("Fp_11"), col("Fp_12")], 1), np.stack([col("Fp_21"), col("Fp_22")], 1)], 1)
        G = F @ np.linalg.inv(Fp)
        phie = np.degrees(np.arctan2(G[:, 1, 0] - G[:, 0, 1], G[:, 0, 0] + G[:, 1, 1]))
        dthp = np.asarray(m.point_data["dtheta_deg"], float).reshape(len(p), -1)[:, 0]
        dth = dthp[tri].mean(1)
        beta = (phie - dth + 180.0) % 360.0 - 180.0
    return dict(area=area, tau=tau, tauc=tauc, phi=phi, xc=xc, T=Tc, press=press, beta=beta, dth=dth)


def stats_for(s, mask, thr):
    """one row of statistics for the elements in mask"""
    w = s["area"][mask]
    tau = s["tau"][:, mask]
    tc = s["tauc"][mask]
    n = int(mask.sum())
    row = {"n_elem": n, "area_um2": float(w.sum())}
    if n == 0:
        return row
    tm = np.abs(tau).max(0)
    q = wquantile(tm, w, [0.5, 0.9, 0.95, 0.99])
    row.update(taumax_max=float(tm.max()), taumax_p99=q[3], taumax_p95=q[2], taumax_p90=q[1], taumax_median=q[0],
               taumax_mean=wmean(tm, w))
    for a in range(3):
        ta = np.abs(tau[a])
        row["tau%d_absmax" % (a + 1)] = float(ta.max())
        row["tau%d_absp99" % (a + 1)] = wquantile(ta, w, [0.99])[0]
        row["tau%d_absmean" % (a + 1)] = wmean(ta, w)
        row["tau%d_min" % (a + 1)] = float(tau[a].min())
        row["tau%d_max" % (a + 1)] = float(tau[a].max())
    qc = wquantile(tc, w, [0.5])
    row.update(tauc_min=float(tc.min()), tauc_median=qc[0], tauc_mean=wmean(tc, w), tauc_max=float(tc.max()))
    r = tm / np.maximum(tc, 1e-12)
    qr = wquantile(r, w, [0.5, 0.99])
    row.update(ratio_median=qr[0], ratio_p99=qr[1], area_frac_tau_ge_tauc=float(w[r >= 1.0].sum() / w.sum()))
    for t in thr:
        row["area_frac_taumax_ge_%g" % t] = float(w[tm >= t].sum() / w.sum())
    if s.get("T") is not None:
        Tm = s["T"][mask]
        qt = wquantile(Tm, w, [0.5, 0.99])
        row.update(T_max=float(Tm.max()), T_p99=qt[1], T_median=qt[0], T_mean=wmean(Tm, w))
    if s.get("press") is not None:
        pr = s["press"][mask]
        qp = wquantile(pr, w, [0.01, 0.5, 0.99])
        row.update(pressure_max=float(pr.max()), pressure_p99=qp[2], pressure_median=qp[1], pressure_p01=qp[0], pressure_min=float(pr.min()))
    if s.get("beta") is not None:
        ab = np.abs(s["beta"][mask])
        qb = wquantile(ab, w, [0.5, 0.99])
        ad = np.abs(s["dth"][mask])
        row.update(beta_absmax_deg=float(ab.max()), beta_absp99_deg=qb[1], beta_absmedian_deg=qb[0],
                   dtheta_absmax_deg=float(ad.max()), dtheta_absp99_deg=wquantile(ad, w, [0.99])[0])
    return row


def selections(s, phisplit, xrange):
    base = np.ones(len(s["area"]), bool)
    if xrange is not None:
        base = (s["xc"] >= xrange[0]) & (s["xc"] <= xrange[1])
    return {"all": base, "bulk": base & (s["phi"] >= phisplit), "core": base & (s["phi"] < phisplit)}


# ------------------------------------------------------------------ one run
def process_run(files, directory, args):
    pvd = read_pvd(directory)
    snaps = []
    for f in files:
        t, tag = snapshot_time(f, pvd, args.dt, args.all_tags)
        if tag is None or (t is None and not args.all_tags):
            continue
        if t is not None and args.tmax is not None and t > args.tmax:
            continue
        snaps.append((t if t is not None else np.nan, tag, f))
    snaps.sort(key=lambda r: (np.isnan(r[0]), r[0]))
    if args.every > 1:
        snaps = snaps[::args.every]
    rows = []
    for t, tag, f in snaps:
        try:
            s = load_snapshot(f)
        except Exception as e:
            sys.stderr.write("skip %s: %s\n" % (os.path.basename(f), e))
            continue
        for name, mask in selections(s, args.phisplit, args.xrange).items():
            r = {"t_ns": t, "snapshot": tag, "selection": name}
            r.update(stats_for(s, mask, args.thr))
            rows.append(r)
        sys.stderr.write("\r  %s  %d/%d   " % (tag, len(rows) // 3, len(snaps)))
    sys.stderr.write("\n")
    return rows


def write_csv(path, rows):
    if not rows:
        return
    cols = []
    for r in rows:
        for k in r:
            if k not in cols:
                cols.append(k)
    with open(path, "w") as fh:
        wr = csv.writer(fh)
        wr.writerow(cols)
        for r in rows:
            wr.writerow([("%.6g" % r[k]) if isinstance(r.get(k), float) else r.get(k, "") for k in cols])


def summarise(tag, rows, window):
    """peaks over time and window means of the 'all' rows"""
    a = [r for r in rows if r["selection"] == "all" and r.get("n_elem")]
    if not a:
        return None
    t = np.array([r["t_ns"] for r in a])
    v, ang = run_label(tag)
    out = {"run": tag, "speed_m_s": v if v is not None else "", "angle_deg": ang if ang is not None else "", "n_snapshots": len(a)}
    for key in ("taumax_max", "taumax_p99", "taumax_median", "tauc_median"):
        x = np.array([r[key] for r in a])
        i = int(np.nanargmax(x))
        out["peak_" + key] = float(x[i])
        out["t_of_peak_" + key] = float(t[i])
    if "T_max" in a[0]:
        Tmx = np.array([r["T_max"] for r in a])
        out["peak_T_max"] = float(Tmx.max())
        out["t_of_peak_T_max"] = float(t[int(np.argmax(Tmx))])
        jump = heating_jump(t, np.array([r["T_mean"] for r in a]))
        if jump:
            out["t_of_largest_mean_T_rise"] = jump[0]
            out["largest_mean_T_rise_K_per_ns"] = jump[1]
        jump = heating_jump(t, Tmx)
        if jump:
            out["t_of_largest_T_max_rise"] = jump[0]
            out["largest_T_max_rise_K_per_ns"] = jump[1]
    if "beta_absmax_deg" in a[0]:
        bm = np.array([r["beta_absmax_deg"] for r in a])
        out["peak_beta_absmax_deg"] = float(bm.max())
        out["t_of_peak_beta_absmax_deg"] = float(t[int(np.argmax(bm))])
        out["peak_beta_absp99_deg"] = float(max(r["beta_absp99_deg"] for r in a))
        dm = np.array([r["dtheta_absmax_deg"] for r in a])
        out["peak_dtheta_absmax_deg_vtu"] = float(dm.max())
        out["t_of_peak_dtheta_absmax_deg_vtu"] = float(t[int(np.argmax(dm))])
    if "pressure_max" in a[0]:
        pm = np.array([r["pressure_p99"] for r in a])
        out["peak_pressure_p99"] = float(pm.max())
        out["t_of_peak_pressure_p99"] = float(t[int(np.argmax(pm))])
    if window is not None:
        sel = (t >= window[0]) & (t <= window[1])
        if sel.any():
            for key in ("taumax_max", "taumax_p99", "taumax_median", "tauc_median"):
                out["window_mean_" + key] = float(np.mean([r[key] for r, k in zip(a, sel) if k]))
            out["window"] = "%g-%g ns (%d snapshots)" % (window[0], window[1], int(sel.sum()))
        else:
            out["window"] = "no snapshot in %g-%g ns" % tuple(window)
    return out


# ------------------------------------------------------------------ figures
def make_plots(allrows, outdir, thr_ref, window, marks, logs=None):
    logs = logs or {}
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        sys.stderr.write("matplotlib missing: no figures\n")
        return
    runs = list(allrows.keys())
    # per run figure(s): tau_max statistics and tau_c, and the three systems
    fig, axs = plt.subplots(len(runs), 4, figsize=(22, 3.6 * len(runs)), squeeze=False)
    for i, tag in enumerate(runs):
        a = [r for r in allrows[tag] if r["selection"] == "all" and r.get("n_elem") and not np.isnan(r["t_ns"])]
        if not a:
            continue
        t = np.array([r["t_ns"] for r in a])
        ax = axs[i, 0]
        ax.plot(t, [r["taumax_max"] for r in a], color="0.2", lw=1.0, label="max")
        ax.plot(t, [r["taumax_p99"] for r in a], color="C3", lw=1.4, label="99th percentile")
        ax.plot(t, [r["taumax_median"] for r in a], color="C0", lw=1.4, label="median")
        ax.plot(t, [r["tauc_median"] for r in a], color="C2", lw=1.4, ls="--", label=r"$\tau_c$ median")
        ax.fill_between(t, [r["tauc_min"] for r in a], [r["tauc_max"] for r in a], color="C2", alpha=0.12, lw=0, label=r"$\tau_c$ range")
        if thr_ref:
            ax.axhline(thr_ref, color="k", ls=":", lw=1, label="%g GPa reference" % thr_ref)
        if window:
            ax.axvspan(window[0], window[1], color="0.9", zorder=0)
        ax.set_ylabel(r"$\max_a|\tau_a|$ per element (GPa)")
        ax.set_xlabel("t (ns)")
        ax.set_title(tag, fontsize=8)
        ax.grid(alpha=0.25)
        if i == 0:
            ax.legend(fontsize=7, ncol=2)
        ax = axs[i, 1]
        for k in (1, 2, 3):
            ax.plot(t, [r["tau%d_absp99" % k] for r in a], lw=1.2, label=r"$\tau_%d$ p99 of $|\tau_%d|$" % (k, k))
        ax.set_xlabel("t (ns)")
        ax.set_ylabel("GPa")
        ax.grid(alpha=0.25)
        ax.set_title("per slip system", fontsize=8)
        if i == 0:
            ax.legend(fontsize=7)
        ax = axs[i, 2]
        if "T_max" in a[0]:
            ax.plot(t, [r["T_max"] for r in a], color="C3", lw=1.3, label="T max")
            ax.plot(t, [r["T_p99"] for r in a], color="C1", lw=1.2, label="T p99")
            ax.plot(t, [r["T_mean"] for r in a], color="C0", lw=1.2, label="T mean")
            ax.set_ylabel("T (K)")
            jump = heating_jump(t, np.array([r["T_mean"] for r in a]))
            if jump:
                ax.axvline(jump[0], color="C3", ls="--", lw=0.8)
                ax.set_title("largest rise of mean T at %.2f ns (%.1f K/ns)" % jump, fontsize=8)
            if "pressure_p99" in a[0]:
                ax2 = ax.twinx()
                ax2.plot(t, [r["pressure_p99"] for r in a], color="0.4", lw=0.9, ls=":", label="pressure p99")
                ax2.set_ylabel("pressure p99 (GPa)", color="0.4")
            if i == 0:
                ax.legend(fontsize=7, loc="upper left")
        ax.set_xlabel("t (ns)")
        ax.grid(alpha=0.25)
        ax = axs[i, 3]
        if "beta_absmax_deg" in a[0]:
            ax.semilogy(t, np.maximum([r["beta_absmax_deg"] for r in a], 1e-6), color="C3", lw=1.3, label=r"$|\beta|$ max")
            ax.semilogy(t, np.maximum([r["beta_absp99_deg"] for r in a], 1e-6), color="C1", lw=1.2, label=r"$|\beta|$ p99")
            ax.semilogy(t, np.maximum([r["beta_absmedian_deg"] for r in a], 1e-6), color="C0", lw=1.2, label=r"$|\beta|$ median")
            ax.set_ylabel(r"mismatch $\beta=\varphi_e-\theta$ (deg)")
            ax2 = ax.twinx()
            ax2.plot(t, [r["dtheta_absmax_deg"] for r in a], color="0.35", lw=1.0, ls=":", label="max|dtheta| (VTU)")
            if tag in logs:
                ax2.plot(logs[tag]["t_ns"], logs[tag]["max_abs_dtheta_deg"], color="0.1", lw=0.9, label="max_abs_dtheta_deg (log)")
            ax2.set_ylabel("rotation (deg)", color="0.35")
            if i == 0:
                ax.legend(fontsize=7, loc="upper left")
                ax2.legend(fontsize=7, loc="lower right")
        ax.set_xlabel("t (ns)")
        ax.grid(alpha=0.25)
        for tm in marks:
            for aa in axs[i]:
                aa.axvline(tm, color="0.5", ls="-.", lw=0.8)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "tau_vs_time.png"), dpi=140)
    plt.close(fig)
    if len(runs) > 1:
        fig, axs = plt.subplots(1, 2, figsize=(12, 4.2))
        for tag in runs:
            a = [r for r in allrows[tag] if r["selection"] == "all" and r.get("n_elem") and not np.isnan(r["t_ns"])]
            if not a:
                continue
            v, ang = run_label(tag)
            lab = ("%s m/s" % v if v is not None else tag) + (" a%s" % ang if ang is not None else "")
            t = np.array([r["t_ns"] for r in a])
            axs[0].plot(t, [r["taumax_p99"] for r in a], lw=1.3, label=lab)
            axs[1].plot(t, [r["taumax_median"] for r in a], lw=1.3, label=lab)
        for ax, ttl in zip(axs, ("99th percentile of $\\max_a|\\tau_a|$ (GPa)", "median of $\\max_a|\\tau_a|$ (GPa)")):
            if thr_ref:
                ax.axhline(thr_ref, color="k", ls=":", lw=1)
            ax.set_xlabel("t (ns)")
            ax.set_ylabel(ttl)
            ax.grid(alpha=0.25)
        axs[0].legend(fontsize=7)
        fig.tight_layout()
        fig.savefig(os.path.join(outdir, "tau_compare.png"), dpi=140)
        plt.close(fig)


# ------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("paths", nargs="+", help="run directories or .vtu files")
    ap.add_argument("--out", default="tau_post")
    ap.add_argument("--dt", type=float, default=0.01, help="ns per loaded step (to turn _s<n> into time)")
    ap.add_argument("--window", type=float, nargs=2, default=None, metavar=("T0", "T1"), help="loading window (ns) for the summary")
    ap.add_argument("--phisplit", type=float, default=0.5, help="cell-mean phi that separates bulk from boundary core")
    ap.add_argument("--thr", type=float, nargs="*", default=[0.1, 0.2, 0.5], help="stress thresholds (GPa) for area fractions")
    ap.add_argument("--xrange", type=float, nargs=2, default=None, help="only elements with centroid x in [X0,X1] (um)")
    ap.add_argument("--every", type=int, default=1, help="use every N-th snapshot")
    ap.add_argument("--tmax", type=float, default=None)
    ap.add_argument("--all-tags", action="store_true", help="also read snapshots whose time cannot be determined")
    ap.add_argument("--ref", type=float, default=0.5, help="reference line in the figures (GPa), 0 = none")
    ap.add_argument("--tend", type=float, default=30.0, help="expected end time (ns); warn if the series stops earlier (0 = no check)")
    ap.add_argument("--mark", type=float, nargs="*", default=[], help="times (ns) to mark in the figures, e.g. --mark 6")
    ap.add_argument("--log", nargs="*", default=[], metavar="CSV",
                    help="solver log CSV(s), one per run path in the same order ('-' to skip one); their t_ns and max_abs_dtheta_deg are overlaid and summarised")
    ap.add_argument("--no-plots", action="store_true")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    allrows, summaries, logs = {}, [], {}
    for ip, p in enumerate(args.paths):
        if os.path.isdir(p):
            files = sorted(glob.glob(os.path.join(p, "*.vtu")))
            directory = p
        else:
            files = [p]
            directory = os.path.dirname(p) or "."
        if not files:
            sys.stderr.write("no .vtu files in %s\n" % p)
            continue
        tag = run_tag(directory, files)
        sys.stderr.write("run %s: %d files\n" % (tag, len(files)))
        rows = process_run(files, directory, args)
        if not rows:
            sys.stderr.write("  nothing usable\n")
            continue
        allrows[tag] = rows
        if ip < len(args.log) and args.log[ip] != "-":
            lg = read_log(args.log[ip])
            if lg is not None:
                logs[tag] = lg
        coverage_report(tag, rows, args.tend if args.tend > 0 else None)
        write_csv(os.path.join(args.out, "tau_stats_%s.csv" % tag), rows)
        s = summarise(tag, rows, args.window)
        if s:
            if tag in logs:
                lg = logs[tag]
                j = int(np.nanargmax(lg["max_abs_dtheta_deg"]))
                s["log_peak_max_abs_dtheta_deg"] = float(lg["max_abs_dtheta_deg"][j])
                s["log_t_of_peak_max_abs_dtheta_deg"] = float(lg["t_ns"][j])
            summaries.append(s)
    write_csv(os.path.join(args.out, "tau_summary.csv"), summaries)
    for s in summaries:
        sys.stdout.write("%s\n  peak max|tau|=%.3f GPa at %.2f ns; peak p99=%.3f GPa at %.2f ns; peak median=%.3f GPa; tau_c median up to %.3f GPa\n" % (
            s["run"], s["peak_taumax_max"], s["t_of_peak_taumax_max"], s["peak_taumax_p99"], s["t_of_peak_taumax_p99"],
            s["peak_taumax_median"], s["peak_tauc_median"]))
        if "peak_beta_absmax_deg" in s:
            sys.stdout.write("  mismatch beta = phi_e - theta: peak max %.2f deg at %.2f ns (p99 peak %.2f deg); max|dtheta| from the VTUs %.1f deg at %.2f ns\n" % (
                s["peak_beta_absmax_deg"], s["t_of_peak_beta_absmax_deg"], s["peak_beta_absp99_deg"],
                s["peak_dtheta_absmax_deg_vtu"], s["t_of_peak_dtheta_absmax_deg_vtu"]))
        if "log_peak_max_abs_dtheta_deg" in s:
            sys.stdout.write("  log max_abs_dtheta_deg: peak %.1f deg at %.2f ns\n" % (s["log_peak_max_abs_dtheta_deg"], s["log_t_of_peak_max_abs_dtheta_deg"]))
        if "t_of_largest_mean_T_rise" in s:
            sys.stdout.write("  largest rise of mean T: %.1f K/ns at %.2f ns; largest rise of max T: %.1f K/ns at %.2f ns; peak T max %.0f K at %.2f ns\n" % (
                s["largest_mean_T_rise_K_per_ns"], s["t_of_largest_mean_T_rise"], s.get("largest_T_max_rise_K_per_ns", float("nan")),
                s.get("t_of_largest_T_max_rise", float("nan")), s["peak_T_max"], s["t_of_peak_T_max"]))
        if "window" in s and "window_mean_taumax_p99" in s:
            sys.stdout.write("  window %s: mean p99=%.3f, mean median=%.3f, mean max=%.3f GPa\n" % (
                s["window"], s["window_mean_taumax_p99"], s["window_mean_taumax_median"], s["window_mean_taumax_max"]))
    if allrows and not args.no_plots:
        make_plots(allrows, args.out, args.ref if args.ref > 0 else None, args.window, args.mark, logs)
    sys.stdout.write("wrote %s\n" % args.out)


if __name__ == "__main__":
    main()
