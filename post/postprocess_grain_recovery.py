#!/usr/bin/env python3
"""postprocess_grain_recovery.py -- grain formation, recovery and retention diagnostics from saved snapshots.

usage:  python3 post/postprocess_grain_recovery.py --config grain_recovery.toml [--run DIR] [--base TAG] [--out DIR] [--np N]
        python3 post/postprocess_grain_recovery.py --config <out>/resolved_config.json      (exact rerun of an earlier analysis)

Reads the run's VTU snapshots (<base>_relaxed.vtu and every frame listed in <base>.pvd), never writes into the run folder.
Reuses grain_count.py (VTU reader, orientation-plateau segmentation) and the field conventions of DATA_FIELDS.md.

Outputs (in --out): grain_recovery.mp4 (three synchronised maps + histories), snapshots.pdf/png, histories.pdf/png,
recovery_<patch>.pdf/png, relaxation.pdf/png, point_history.csv, patch_history.csv, selection.csv, summary.txt,
resolved_config.json.  See README_grain_recovery.md.
"""
import argparse, json, os, re, sys, math, hashlib, subprocess, shutil, time
from multiprocessing import Pool
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import grain_count as gc   # noqa: E402  (VTU reader and segmentation of the repository)

try:
    import tomllib
except ImportError:  # python < 3.11
    tomllib = None

# ----------------------------------------------------------------------------------------------------------------- config
DEFAULTS = {
    "run": {"dir": "", "base": "", "tmin": None, "tmax": None, "nominal_dt_ns": 0.1, "include_relaxed": True},
    "fields": {"theta": "theta_lattice_deg", "phi": "phi", "rhoS": "rhoS", "rhoG": "rhoG", "slip_path": "slip_path_phys", "slip_rate_path": "slip_path_eff",
               "label": "Label", "u": "u", "sigma": ["sig_11", "sig_22", "sig_12", "sig_21", "sig_33"],
               # optional saved rate fields (none are written by the 7.x solver; give names here if a later version saves them)
               "rate_storage": "", "rate_dynrec": "", "rate_sweep": ""},
    "segmentation": {"grad_deg_per_um": 2.0, "merge_deg": 2.0, "minarea_um2": 0.1, "sym_deg": 90.0},
    "selection": {"mode": "centers", "center_method": "inscribed", "min_grain_area_um2": 0.1, "max_grains": 60, "n_new_grains": 4, "new_grain_dtheta_deg": 5.0, "min_interior_dist_um": 0.15,
                  "near_boundary_dist_um": 0.05, "patch_radius_um": 0.2, "unrefined": True, "points": []},
    "model": {"K": 0.05, "d_um": 5e-4, "b_um": 2.556e-4, "C_D": None, "C_A_um": None, "recov": None, "rho_floor": 1e-3,
              "rho0": 10.0},
    "scales": {"theta_deg": None, "rhoS_um2": None, "rhoG_um2": None, "vis_floor_um2": 1e-2, "theta_pad_deg": 20.0},
    "timeline": {"markers_ns": [3.0, 13.0], "marker_labels": ["loading ends", "~99% of W converted to heat"],
                 "snapshot_times_ns": [0.0, 3.0, 13.0, 30.0]},
    "render": {"fps": 10, "dpi": 110, "phi_panel": False, "deformed": False, "mp4": True},
    "output": {"dir": "grain_recovery_out"},
}


def deep_update(a, b):
    for k, v in b.items():
        if isinstance(v, dict) and isinstance(a.get(k), dict):
            deep_update(a[k], v)
        else:
            a[k] = v
    return a


def load_config(path):
    cfg = json.loads(json.dumps(DEFAULTS))
    if path:
        if path.endswith(".json"):
            user = json.load(open(path))
        else:
            if tomllib is None:
                sys.exit("TOML needs python >= 3.11 (or give a .json config)")
            user = tomllib.load(open(path, "rb"))
        deep_update(cfg, user)
    return cfg


# ----------------------------------------------------------------------------------------------------------------- frames
def frame_list(run, base, include_relaxed, notes):
    """(t_ns, path) of every snapshot: the relaxed frame (t = 0) and every frame of <base>.pvd that exists on disk.
    Missing files, duplicate times and gaps against the nominal spacing are recorded in notes."""
    fs = []
    rel = os.path.join(run, base + "_relaxed.vtu")
    if include_relaxed and os.path.exists(rel):
        fs.append((0.0, rel, "relaxed"))
    pvd = os.path.join(run, base + ".pvd")
    if not os.path.exists(pvd):
        sys.exit(f"no {pvd}")
    listed = re.findall(r'timestep="([^"]+)"[^>]*file="([^"]+)"', open(pvd).read())
    notes["pvd_entries"] = len(listed)
    missing = []
    for t, f in listed:
        p = os.path.join(run, os.path.basename(f))
        if "relaxed" in p:
            continue
        if os.path.exists(p):
            fs.append((float(t), p, "pvd"))
        else:
            missing.append((float(t), os.path.basename(f)))
    seen = set(); uniq = []
    for t, f in missing:
        if f not in seen: uniq.append((t, f)); seen.add(f)
    notes["listed_but_missing"] = uniq
    fs.sort(key=lambda a: a[0])
    out, dup = [], []
    for t, p, src in fs:
        if out and abs(t - out[-1][0]) < 1e-9:
            dup.append((t, os.path.basename(p), os.path.basename(out[-1][1])))
            continue
        out.append((t, p, src))
    notes["duplicates_dropped"] = dup
    return out


def nominal_gaps(times, dt, notes):
    if len(times) < 2:
        return
    exp = np.arange(round(times[0] / dt) * dt, times[-1] + 0.5 * dt, dt)
    have = np.array(times)
    gaps = [float(round(e, 6)) for e in exp if np.min(np.abs(have - e)) > 0.25 * dt]
    notes["nominal_times_absent"] = gaps
    off = [float(t) for t in have if abs(t / dt - round(t / dt)) > 0.25 and t > 0]
    notes["off_grid_times"] = off   # e.g. a final frame at t_end or a TSTOP frame


# ----------------------------------------------------------------------------------------------------------------- mesh
class Geo:
    """Fixed reference mesh: triangles, areas, centroids, edges, point location and P1 gradients."""

    def __init__(self, pts, tri, th0, phi0):
        self.P, self.T = pts, tri
        self.M = gc.Mesh(pts, tri, th0, phi0)           # edges, graph, areas (grain_count)
        self.A, self.cen = self.M.A, self.M.cen
        p0, p1, p2 = pts[tri[:, 0]], pts[tri[:, 1]], pts[tri[:, 2]]
        self.e1, self.e2 = p1 - p0, p2 - p0
        self.det = self.e1[:, 0] * self.e2[:, 1] - self.e1[:, 1] * self.e2[:, 0]
        import matplotlib.tri as mtri
        self.tr = mtri.Triangulation(pts[:, 0], pts[:, 1], tri)
        self.finder = self.tr.get_trifinder()
        self.h = self.M.h
        ed = np.sort(np.concatenate([tri[:, [0, 1]], tri[:, [1, 2]], tri[:, [2, 0]]]), axis=1)
        u, c = np.unique(ed, axis=0, return_counts=True)
        self.outer = pts[u[c == 1]]                               # domain boundary edges (one triangle), (n, 2, 2)
        self.sig = hashlib.sha1(np.ascontiguousarray(pts.astype(np.float32)).tobytes() + np.ascontiguousarray(tri.astype(np.int64)).tobytes()).hexdigest()

    def grad(self, f):
        """element gradient of a P1 field (per triangle, units of f per um)."""
        t = self.T
        d1, d2 = f[t[:, 1]] - f[t[:, 0]], f[t[:, 2]] - f[t[:, 0]]
        gx = (d1 * self.e2[:, 1] - d2 * self.e1[:, 1]) / self.det
        gy = (-d1 * self.e2[:, 0] + d2 * self.e1[:, 0]) / self.det
        return gx, gy

    def locate(self, x, y):
        k = int(self.finder(x, y))
        if k < 0:
            raise ValueError(f"point ({x}, {y}) is outside the mesh")
        p0 = self.P[self.T[k, 0]]
        r = np.array([x, y]) - p0
        l1 = (r[0] * self.e2[k, 1] - r[1] * self.e2[k, 0]) / self.det[k]
        l2 = (-r[0] * self.e1[k, 1] + r[1] * self.e1[k, 0]) / self.det[k]
        return k, np.array([1 - l1 - l2, l1, l2])

    def segment(self, th, phi, par):
        M = self.M
        M.th, M.phi = th, phi
        d1, d2 = th[self.T[:, 1]] - th[self.T[:, 0]], th[self.T[:, 2]] - th[self.T[:, 0]]
        M.g = np.hypot((d1 * self.e2[:, 1] - d2 * self.e1[:, 1]) / self.det, (-d1 * self.e2[:, 0] + d2 * self.e1[:, 0]) / self.det)
        M.tm = th[self.T].mean(1)
        return gc.segment(M, par["grad_deg_per_um"], par["minarea_um2"], par["merge_deg"])

    def boundary_edges(self, terr):
        M = self.M
        ta, tb = terr[M.pa], terr[M.pb]
        e = np.flatnonzero((ta != tb) & (ta >= 0) & (tb >= 0))
        return self.P[M.en[e]]                               # (n, 2, 2) segments in reference coordinates


def seg_dist(p, segs):
    """distance from point p to the nearest of the segments (n, 2, 2)."""
    if len(segs) == 0:
        return np.inf
    a, b = segs[:, 0], segs[:, 1]
    ab = b - a
    t = np.clip(((p - a) * ab).sum(1) / np.maximum((ab * ab).sum(1), 1e-30), 0, 1)
    q = a + t[:, None] * ab
    return float(np.sqrt(((q - p) ** 2).sum(1)).min())


def seg_dist_many(P, segs):
    from scipy.spatial import cKDTree
    if len(segs) == 0:
        return np.full(len(P), np.inf)
    mids = 0.5 * (segs[:, 0] + segs[:, 1])
    L = np.sqrt(((segs[:, 1] - segs[:, 0]) ** 2).sum(1)).max()
    d0, _ = cKDTree(mids).query(P)
    return np.maximum(d0 - 0.5 * L, 0.0)                     # lower bound within half an edge length; exact enough for selection


# ----------------------------------------------------------------------------------------------------------------- reading
def read_frame(path, F):
    d = gc.read_vtu(path)
    pts = d["points"][:, :2].astype(float)
    conn, off, typ = d["connectivity"], d["offsets"], d["types"]
    start = np.r_[0, off[:-1]]
    istri = typ == 5
    tri = np.array([conn[s:s + 3] for s, t in zip(start, typ) if t == 5])

    def cell(name):
        return None if not name or name not in d else np.asarray(d[name], float).ravel()[istri] if len(np.asarray(d[name]).ravel()) == len(typ) else np.asarray(d[name], float).ravel()

    def point(name):
        return None if not name or name not in d else np.asarray(d[name], float).ravel()

    out = dict(pts=pts, tri=tri, theta=point(F["theta"]), phi=point(F["phi"]), rhoS=cell(F["rhoS"]), rhoG=cell(F["rhoG"]),
               slip=cell(F["slip_path"]), slipE=cell(F.get("slip_rate_path", "")), label=cell(F["label"]))
    out["u"] = d[F["u"]][:, :2].astype(float) if F.get("u") and F["u"] in d else None
    sg = [cell(s) for s in F["sigma"]]
    out["sigma"] = None if any(s is None for s in sg) else np.stack(sg, 1)
    for k in ("rate_storage", "rate_dynrec", "rate_sweep"):
        out[k] = cell(F.get(k, ""))
    for k in ("theta", "phi", "rhoS", "rhoG"):
        if out[k] is None:
            raise KeyError(f"{os.path.basename(path)}: field '{F[k]}' ({k}) not found")
    return out


# ----------------------------------------------------------------------------------------------------------------- run parameters
def run_parameters(run, base, cfg, notes):
    """SSD-law parameters of the run that produced the data: <base>_meta.json, then <base>.args, then the option dump in
    <base>.out; constants that are not options (K, d, b, floor) come from the configuration (documented values of the
    implementation document, section 5.4) and are flagged as such."""
    m = dict(cfg["model"])
    src = {}
    opts = {}
    meta = os.path.join(run, base + "_meta.json")
    if os.path.exists(meta):
        try:
            j = json.load(open(meta))
            flat = {}

            def walk(o, pre=""):
                if isinstance(o, dict):
                    for k, v in o.items():
                        walk(v, k)
                else:
                    flat[pre] = o
            walk(j)
            for k in ("cd", "ca", "recov", "C_D", "C_A"):
                if k in flat:
                    opts[k] = flat[k]; src[k] = "meta.json"
        except Exception as e:  # noqa: BLE001
            notes.setdefault("warnings", []).append(f"meta.json unreadable: {e}")
    for fn, tag in ((os.path.join(run, base + ".args"), ".args"), (os.path.join(run, base + ".out"), ".out")):
        if not os.path.exists(fn):
            continue
        txt = open(fn, errors="ignore").read()
        for k in ("cd", "ca", "recov"):
            if k in opts:
                continue
            mm = re.search(r"-%s\s+([0-9.eE+-]+)" % k, txt) or re.search(r"name:-%s value: ([0-9.eE+-]+)" % k, txt)
            if mm:
                opts[k] = float(mm.group(1)); src[k] = tag
        mv = re.search(r"code version\s*:?\s*([0-9.]+)", txt)
        if mv:
            notes["code_version"] = mv.group(1)
    if m["C_D"] is None:
        m["C_D"] = float(opts.get("cd", opts.get("C_D", 100.0))); src.setdefault("cd", "default (not found in the run files)")
    if m["C_A_um"] is None:
        m["C_A_um"] = float(opts.get("ca", opts.get("C_A", math.sqrt(10)))); src.setdefault("ca", "default (not found in the run files)")
    if m["recov"] is None:
        m["recov"] = int(opts.get("recov", 2)); src.setdefault("recov", "default 2 (not found in the run files)")
    src["K, d, b, floor, rho0"] = "configuration [model] (documented constants, cosserat_kwc_production.md section 5.4; not options of the solver)"
    notes["model_parameter_sources"] = src
    if "code_version" not in notes:
        notes["code_version"] = "not recorded in the run files (runs before version 7.1)"
    return m


# ----------------------------------------------------------------------------------------------------------------- selection
def select(geo, F0, FL, segL, cfg, notes):
    """Points and patches, chosen on the final frame (FL) and fixed in reference coordinates."""
    S = cfg["selection"]
    sym = cfg["segmentation"]["sym_deg"]
    rows, bnd = gc.analyse(geo.M, segL, (FL["label"] if FL["label"] is not None else np.zeros(len(geo.T))).astype(int), sym)
    terr, plat, mu = segL["terr"], segL["plat"], segL["mu"]
    tm0 = F0["theta"][geo.T].mean(1)
    tmL = FL["theta"][geo.T].mean(1)
    segs = geo.boundary_edges(terr)
    # depth used for choosing interior points: distance to the nearest grain boundary OR domain edge
    dist = np.minimum(seg_dist_many(geo.cen, segs), seg_dist_many(geo.cen, geo.outer))
    grains = []
    for r in rows:
        g = r["id"]
        p = plat == g
        a = geo.A[p]
        d_init = float((a * (tmL[p] - tm0[p])).sum() / a.sum())          # rotation of this material since the first frame
        grains.append(dict(id=g, area=r["plateau_area"], theta=r["theta_deg"], orig_label=r["orig_label"], d_init=d_init,
                           maxdist=float(dist[p].max()) if p.any() else 0.0))
    # parent of an initial grain = the final grain of that label that kept its initial orientation best (largest area on ties):
    # the remnant of the original grain, whatever its size
    for g in grains:
        same = [q for q in grains if q["orig_label"] == g["orig_label"]]
        parent = min(same, key=lambda q: (round(abs(q["d_init"]), 1), -q["area"]))
        g["is_parent"] = parent["id"] == g["id"]
        g["d_parent"] = 0.0 if g["is_parent"] else float(abs(g["theta"] - parent["theta"]))
    thr = S["new_grain_dtheta_deg"]
    new = [g for g in grains if (not g["is_parent"] and g["d_parent"] >= thr) or abs(g["d_init"]) >= thr]
    new.sort(key=lambda g: -g["area"])
    notes["final_grains"] = len(grains)
    notes["new_grain_candidates"] = len(new)
    sel = []

    def add(kind, label, k, rationale, grain=None, radius=None):
        x, y = geo.cen[k]
        R = radius if radius is not None else min(S["patch_radius_um"], max(0.5 * dist[k], 1.01 * geo.h))
        mem = np.flatnonzero(np.hypot(*(geo.cen - geo.cen[k]).T) <= R)
        if len(mem) < 3:
            mem = np.argsort(np.hypot(*(geo.cen - geo.cen[k]).T))[:3]
        sel.append(dict(id=f"P{len(sel) + 1}", kind=kind, label=label, x=float(x), y=float(y), patch_radius_um=float(R),
                        patch_tris=mem.tolist(), patch_area_um2=float(geo.A[mem].sum()), final_boundary_dist_um=float(dist[k]),
                        final_grain=grain if grain is not None else int(terr[k]), rationale=rationale))

    if S["mode"] == "centers":
        # one point per final grain at its centre: "inscribed" = the plateau triangle farthest from the final boundaries
        # (centre of the largest inscribed circle, always inside the grain); "centroid" = the area-weighted centroid of the
        # plateau, replaced by the inscribed centre when it does not lie inside the grain
        cand = sorted([g for g in grains if g["area"] >= S["min_grain_area_um2"]], key=lambda g: -g["area"])
        if len(cand) > S["max_grains"]:
            notes.setdefault("warnings", []).append(f"{len(cand)} grains above min_grain_area_um2; only the {S['max_grains']} largest are tracked")
            cand = cand[: S["max_grains"]]
        for g in cand:
            p = np.flatnonzero(plat == g["id"])
            k = p[np.argmax(dist[p])]
            how = "deepest point (largest inscribed circle)"
            if S["center_method"] == "centroid":
                c = (geo.cen[p] * geo.A[p, None]).sum(0) / geo.A[p].sum()
                kc = int(geo.finder(c[0], c[1]))
                if kc >= 0 and terr[kc] == g["id"]:
                    k, how = kc, "area-weighted centroid of the plateau"
                else:
                    how = "deepest point (the centroid lies outside the grain)"
            new_ = (not g["is_parent"] and g["d_parent"] >= thr) or abs(g["d_init"]) >= thr
            kind = "center-new" if new_ else "center-kept"
            why = (f"final grain {g['id']} (plateau {g['area']:.3f} um^2, theta {g['theta']:.1f} deg, initial grain {g['orig_label']}): "
                   f"{how}; rotated {g['d_init']:+.1f} deg since the first frame, {g['d_parent']:.1f} deg from the remnant of its initial grain"
                   f" -> {'new' if new_ else 'kept its orientation'} (threshold {thr} deg, label only)")
            if dist[k] < S["min_interior_dist_um"]:
                why += f" [WARNING: centre only {dist[k]:.3f} um from a boundary]"
            add(kind, f"grain {g['id']}" + (" (new)" if new_ else ""), k, why, g["id"])
            sel[-1]["id"] = f"G{g['id']}"
    if S["mode"] in ("auto", "both"):
        for g in new[: S["n_new_grains"]]:
            p = np.flatnonzero(plat == g["id"])
            k = p[np.argmax(dist[p])]
            why = (f"final grain {g['id']} (plateau {g['area']:.3f} um^2, theta {g['theta']:.1f} deg) inside initial grain "
                   f"{g['orig_label']}: rotated {g['d_init']:+.1f} deg since the first frame, {g['d_parent']:.1f} deg from the remnant of "
                   f"that initial grain (the final grain of the same label that kept its orientation); point = plateau triangle "
                   f"farthest from the final boundaries")
            if dist[k] < S["min_interior_dist_um"]:
                why += f" [WARNING: only {dist[k]:.3f} um from a boundary, below min_interior_dist_um]"
            add("new-interior", f"new grain {g['id']}", k, why, g["id"])
            t = np.flatnonzero(terr == g["id"])
            kb = t[np.argmin(np.abs(dist[t] - S["near_boundary_dist_um"]))]
            add("new-near-boundary", f"new grain {g['id']} edge", kb,
                f"same final grain, {dist[kb]:.3f} um from its boundary (target {S['near_boundary_dist_um']} um)", g["id"],
                radius=max(min(S["patch_radius_um"], 0.5 * S["near_boundary_dist_um"]), 1.01 * geo.h))
        if S["unrefined"]:
            old = [g for g in grains if g["is_parent"] and abs(g["d_init"]) < thr]
            if old:
                g = max(old, key=lambda g: g["maxdist"])
                p = np.flatnonzero(plat == g["id"])
                k = p[np.argmax(dist[p])]
                add("unrefined", f"parent grain {g['id']}", k,
                    f"largest-interior grain that kept its initial orientation (rotated {g['d_init']:+.1f} deg); comparison point", g["id"])
            else:
                notes.setdefault("warnings", []).append("no unrefined grain found (every final grain rotated by more than the threshold)")
    if S["mode"] in ("manual", "both"):
        for q in S["points"]:
            k, _ = geo.locate(q["x"], q["y"])
            add(q.get("kind", "manual"), q.get("label", "manual"), k, q.get("rationale", "manual selection"),
                radius=q.get("patch_radius_um"))
            sel[-1]["x"], sel[-1]["y"] = float(q["x"]), float(q["y"])    # keep the given coordinate (not the centroid)
    if not sel:
        sys.exit("no points selected: check [selection] (no new grain found? set mode = 'manual' and give points)")
    for s in sel:
        s["tri"], s["bary"] = geo.locate(s["x"], s["y"])
        s["bary"] = s["bary"].tolist()
    notes["grains_table"] = grains
    return sel


# ----------------------------------------------------------------------------------------------------------------- per-frame pass
G = {}   # worker globals


def init_worker(state):
    G.update(state)


def eval_frame(args):
    i, t, path = args
    geo, cfg, sel = G["geo"], G["cfg"], G["sel"]
    F = read_frame(path, cfg["fields"])
    chk = hashlib.sha1(np.ascontiguousarray(F["pts"].astype(np.float32)).tobytes() + np.ascontiguousarray(F["tri"].astype(np.int64)).tobytes()).hexdigest()
    if chk != geo.sig:
        return dict(i=i, t=t, error="mesh differs from the first frame: no material correspondence by node/element index")
    seg = geo.segment(F["theta"], F["phi"], cfg["segmentation"])
    terr = seg["terr"] if seg else np.full(len(geo.T), -1)
    segs = geo.boundary_edges(terr)
    gx, gy = geo.grad(np.radians(F["theta"]))
    gth2 = gx * gx + gy * gy
    phie = F["phi"][geo.T].mean(1)
    the = F["theta"][geo.T].mean(1)
    pts = []
    for s in sel:
        k, w = s["tri"], np.array(s["bary"])
        nodes = geo.T[k]
        pts.append(dict(theta=float(w @ F["theta"][nodes]), phi=float(w @ F["phi"][nodes]), rhoS=float(F["rhoS"][k]),
                        rhoG=float(F["rhoG"][k]), rhoF=float(F["rhoS"][k] + F["rhoG"][k]), slip=float(F["slip"][k]) if F["slip"] is not None else np.nan,
                        slipE=float(F["slipE"][k]) if F["slipE"] is not None else np.nan,
                        bdist=seg_dist(np.array([s["x"], s["y"]]), segs), grain=int(terr[k]),
                        grain_theta=float(seg["mu"][terr[k]]) if seg and terr[k] >= 0 else np.nan))
    patches = []
    for s in sel:
        m = np.array(s["patch_tris"]); a = geo.A[m]; W = a.sum()
        row = {}
        for nm, v in (("theta", the[m]), ("phi", phie[m]), ("rhoS", F["rhoS"][m]), ("rhoG", F["rhoG"][m]),
                      ("rhoF", F["rhoS"][m] + F["rhoG"][m]),
                      ("slipE", F["slipE"][m] if F["slipE"] is not None else np.full(len(m), np.nan))):
            mu_ = float((a * v).sum() / W); row[nm] = mu_; row[nm + "_sd"] = float(np.sqrt((a * (v - mu_) ** 2).sum() / W))
        # per-triangle quantities kept for the rate reconstruction (patch only)
        row["_rho"] = F["rhoS"][m].astype(float); row["_slip"] = (F["slip"][m] if F["slip"] is not None else np.zeros(len(m))).astype(float)
        row["_phi"] = phie[m].astype(float); row["_gth2"] = gth2[m].astype(float)
        for k in ("rate_storage", "rate_dynrec", "rate_sweep"):
            if F[k] is not None:
                row["_" + k] = F[k][m].astype(float)
        patches.append(row)
    st = {}
    if F["sigma"] is not None:
        s_ = F["sigma"]
        fro = np.sqrt((s_ ** 2).sum(1))
        st = dict(sig_rms=float(np.sqrt((geo.A * fro ** 2).sum() / geo.A.sum())), sig_p99=float(np.percentile(fro, 99)))
    lim = {}
    for k in ("rhoS", "rhoG"):
        v = F[k]
        lim[k] = dict(min=float(v.min()), max=float(v.max()), pos_min=float(v[v > 0].min()) if (v > 0).any() else np.nan,
                      n_neg=int((v < 0).sum()), n_zero=int((v == 0).sum()))
    lim["theta"] = dict(min=float(F["theta"].min()), max=float(F["theta"].max()))
    return dict(i=i, t=t, path=path, slipE_all=(F["slipE"].astype(np.float32) if F["slipE"] is not None else None), n_grains=int(seg["G"]) if seg else 0, segs=segs.astype(np.float32), points=pts, patches=patches,
                stress=st, lim=lim)


# ----------------------------------------------------------------------------------------------------------------- rates
def reconstruct_rates(R, sel, model):
    """Snapshot estimates of storage, dynamic recovery and sweeping recovery per patch (area-weighted), between saved frames."""
    K, d, b, CD, CA, mode, floor = model["K"], model["d_um"], model["b_um"], model["C_D"], model["C_A_um"], model["recov"], model["rho_floor"]
    out = []
    for j, s in enumerate(sel):
        a = np.asarray(G_geoA[s["id"]])
        cum = dict(store=0.0, dyn=0.0, sweep=0.0, pred=0.0, obs=0.0)
        rows = []
        for n in range(len(R)):
            P = R[n]["patches"][j]
            if n == 0:
                rows.append(dict(t=R[n]["t"], **{k: 0.0 for k in cum}, exact=False)); continue
            Q = R[n - 1]["patches"][j]
            r0, r1 = Q["_rho"], P["_rho"]
            dG = np.maximum(P["_slip"] - Q["_slip"], 0.0)
            store = K / b * 0.5 * (np.sqrt(np.maximum(r0, 0)) + np.sqrt(np.maximum(r1, 0))) * dG
            dyn = 2 * d / b * 0.5 * (r0 + r1) * dG
            dphi = P["_phi"] - Q["_phi"]
            if mode == 0:
                sw = CD * np.abs(dphi) * r0
            elif mode == 1:
                sw = CD * np.maximum(dphi, 0) * r0
            else:
                sw = CD * np.maximum(dphi, 0) * np.tanh(CA ** 2 * 0.5 * (P["_gth2"] + Q["_gth2"])) * r0
            sw = np.minimum(sw, np.maximum(r0 + store - dyn - floor, 0))   # the update cannot go below the floor
            exact = all(("_" + k) in P for k in ("rate_storage", "rate_dynrec", "rate_sweep"))
            dt = R[n]["t"] - R[n - 1]["t"]
            if exact:   # saved rates (per ns) of the end frame times the interval: still a rectangle rule, but per-element exact rates
                store, dyn, sw = P["_rate_storage"] * dt, P["_rate_dynrec"] * dt, P["_rate_sweep"] * dt
            W = a.sum()
            inc = dict(store=float((a * store).sum() / W), dyn=float((a * dyn).sum() / W), sweep=float((a * sw).sum() / W),
                       obs=float((a * (r1 - r0)).sum() / W))
            inc["pred"] = inc["store"] - inc["dyn"] - inc["sweep"]
            for k in cum:
                cum[k] += inc[k]
            rows.append(dict(t=R[n]["t"], exact=exact, **{k: cum[k] for k in cum}))
        out.append(rows)
    return out


G_geoA = {}


# ----------------------------------------------------------------------------------------------------------------- rendering
COLORS = ["#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#ff7f0e", "#8c564b", "#e377c2", "#17becf", "#bcbd22", "#7f7f7f"]


def style_for(sel):
    if len(sel) <= len(COLORS):
        return {s["id"]: COLORS[i] for i, s in enumerate(sel)}
    import matplotlib.cm as cm
    cmap = cm.get_cmap("tab20", 20) if len(sel) <= 20 else cm.get_cmap("gist_ncar", len(sel) + 2)
    return {s["id"]: cmap(i) for i, s in enumerate(sel)}


def render_frame(args):
    n, outpng = args
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    from matplotlib.collections import LineCollection
    geo, cfg, sel, R, H, L = G["geo"], G["cfg"], G["sel"], G["R"], G["H"], G["L"]
    col = style_for(sel)
    F = read_frame(R[n]["path"], cfg["fields"])
    t = R[n]["t"]
    panels = [("theta", r"lattice orientation $\theta_t$ (deg)"), ("rhoS", r"SSD density $\rho_S$ ($\mu$m$^{-2}$)"),
              ("rhoG", r"GND density $\rho_G$ ($\mu$m$^{-2}$)"),
              ("slip_rate", r"accumulated slip rate $\Delta(\sum_\alpha\int|\dot\gamma_\alpha|dt)/\Delta t$ (1/ns)")]
    if cfg["render"]["phi_panel"]:
        panels.append(("phi", r"order parameter $\phi$"))
    npn = len(panels)
    fig = plt.figure(figsize=(5.2 * npn, 8.6))
    gsp = fig.add_gridspec(2, 1, height_ratios=[1.55, 1], hspace=0.34)
    top = gsp[0].subgridspec(1, npn, wspace=0.08)
    bot = gsp[1].subgridspec(1, 6, wspace=0.5)
    tr = geo.tr
    if cfg["render"].get("_deformed") and F["u"] is not None:
        import matplotlib.tri as mtri
        tr = mtri.Triangulation(geo.P[:, 0] + F["u"][:, 0], geo.P[:, 1] + F["u"][:, 1], geo.T)
    floor = cfg["scales"]["vis_floor_um2"]
    for c, (key, title) in enumerate(panels):
        ax = fig.add_subplot(top[c])
        if key == "theta":
            pc = ax.tripcolor(tr, F["theta"], shading="gouraud", cmap="twilight", vmin=L["theta"][0], vmax=L["theta"][1])
        elif key == "phi":
            pc = ax.tripcolor(tr, F["phi"], shading="gouraud", cmap="viridis", vmin=0, vmax=1)
        elif key == "slip_rate":
            rf = os.path.join(G["tmpd"], f"rate_{n:05d}.npy")
            v = np.load(rf) if os.path.exists(rf) else np.zeros(len(geo.T))
            v = np.maximum(v, L["slip_rate"][0])                      # display: zero rate drawn at the lower limit
            pc = ax.tripcolor(tr, facecolors=v, cmap="turbo", norm=LogNorm(*L["slip_rate"]))
        else:
            v = np.where(F[key] > floor, F[key], floor)               # display only: raw values unchanged in the CSVs
            pc = ax.tripcolor(tr, facecolors=v, cmap="magma", norm=LogNorm(*L[key]))
        segs = R[n]["segs"]
        if cfg["render"].get("_deformed") and F["u"] is not None and len(segs):
            # map boundary edges through the placement with the nodal displacement of the edge end points
            from scipy.spatial import cKDTree
            idx = G["kd"].query(segs.reshape(-1, 2))[1]
            segs = segs + F["u"][idx].reshape(segs.shape).astype(np.float32)
        ax.add_collection(LineCollection(segs, colors="k", linewidths=0.45, alpha=0.8))
        for s in sel:
            x, y = s["x"], s["y"]
            if cfg["render"].get("_deformed") and F["u"] is not None:
                k, w = s["tri"], np.array(s["bary"]); x, y = x + w @ F["u"][geo.T[k], 0], y + w @ F["u"][geo.T[k], 1]
            ax.plot(x, y, "o", ms=6, mfc=col[s["id"]], mec="w", mew=1.0)
            ax.annotate(s["id"], (x, y), xytext=(4, 4), textcoords="offset points", fontsize=8, color="w",
                        bbox=dict(boxstyle="round,pad=0.1", fc=col[s["id"]], ec="none", alpha=0.85))
        ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
        ax.set_xlim(L["xlim"]); ax.set_ylim(L["ylim"])
        ax.set_title(title, fontsize=10)
        cb = fig.colorbar(pc, ax=ax, orientation="horizontal", fraction=0.05, pad=0.03)
        cb.ax.tick_params(labelsize=8)
        if key in ("rhoS", "rhoG"):
            cb.set_label(f"log scale; values <= {floor:g} drawn at {floor:g} (display floor)", fontsize=7)
        if key == "slip_rate":
            cb.set_label("log scale; average since the previous snapshot" + ("; no rate at the first frame" if n == 0 else ""), fontsize=7)
    geo_lbl = "deformed geometry (x = X + u)" if cfg["render"].get("_deformed") else "reference coordinates X (um)"
    fig.suptitle(f"t = {t:.2f} ns   [{geo_lbl}; boundaries: orientation plateaus, grad < {cfg['segmentation']['grad_deg_per_um']} deg/um]",
                 fontsize=11, y=0.985)
    hist = [("theta", r"$\theta_t$ (deg)", False), ("rhoS", r"$\rho_S$ ($\mu$m$^{-2}$)", True), ("rhoG", r"$\rho_G$ ($\mu$m$^{-2}$)", True),
            ("rhoF", r"$\rho_F$ ($\mu$m$^{-2}$)", True), ("rhoF_pct", r"$\rho_F$ (% of max)", False),
            ("slip_rate", r"slip rate (1/ns; 0 not drawn)", True)]
    htitle = {"rhoF": r"forest $\rho_F=\rho_S+\rho_G$", "rhoF_pct": r"$\rho_F$ / its max over all t",
              "slip_rate": "accumulated slip rate"}
    for c, (key, lab, lg) in enumerate(hist):
        ax = fig.add_subplot(bot[c])
        for s in sel:
            y = H[s["id"]][key]
            ax.plot(H["t"], y, "-", color=col[s["id"]], lw=1.1, label=s["id"])
            ax.plot(H["t"][n], y[n], "o", color=col[s["id"]], ms=3.5)
        if lg:
            ax.set_yscale("log")
        for m, ml in zip(cfg["timeline"]["markers_ns"], cfg["timeline"]["marker_labels"]):
            if H["t"][0] <= m <= H["t"][-1]:
                ax.axvline(m, color="0.4", ls="--", lw=0.8)
        ax.axvline(t, color="k", lw=1.0)
        ax.set_xlabel("t (ns)", fontsize=9); ax.set_ylabel(lab, fontsize=9); ax.tick_params(labelsize=8)
        ax.set_xlim(H["t"][0], H["t"][-1])
        if key in htitle:
            ax.set_title(htitle[key], fontsize=8)
        if key == "rhoF_pct":
            ax.set_ylim(0, 105)
        if c == 0:
            ax.legend(fontsize=7, ncol=2, loc="best", frameon=False)
    mk = ", ".join(f"{m:g} ns: {l}" for m, l in zip(cfg["timeline"]["markers_ns"], cfg["timeline"]["marker_labels"])
                   if H["t"][0] <= m <= H["t"][-1])
    fig.text(0.5, 0.005, ("dashed: " + mk + ";  " if mk else "") + "solid: current frame (saved snapshot, no interpolation)",
             ha="center", fontsize=8)
    fig.savefig(outpng, dpi=cfg["render"]["dpi"])
    plt.close(fig)
    return outpng


# ----------------------------------------------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=None); ap.add_argument("--run"); ap.add_argument("--base"); ap.add_argument("--out")
    ap.add_argument("--np", type=int, default=int(os.environ.get("NP", 1)))
    ap.add_argument("--no-mp4", action="store_true")
    a = ap.parse_args()
    cfg = load_config(a.config)
    if a.run: cfg["run"]["dir"] = a.run
    if a.base: cfg["run"]["base"] = a.base
    if a.out: cfg["output"]["dir"] = a.out
    if a.no_mp4: cfg["render"]["mp4"] = False
    run, base, out = cfg["run"]["dir"], cfg["run"]["base"], cfg["output"]["dir"]
    if not base:
        pv = [f for f in os.listdir(run) if f.endswith(".pvd")]
        if len(pv) != 1:
            sys.exit("give [run] base (several or no .pvd files in the run folder)")
        base = cfg["run"]["base"] = pv[0][:-4]
    if os.path.abspath(out).startswith(os.path.abspath(run) + os.sep) or os.path.abspath(out) == os.path.abspath(run):
        sys.exit("the output folder must not be inside the run folder (the analysis never writes into simulation files)")
    os.makedirs(out, exist_ok=True)
    notes = {"started": time.strftime("%Y-%m-%d %H:%M:%S"), "warnings": []}
    T0 = time.time()

    # frames
    FR = frame_list(run, base, cfg["run"]["include_relaxed"], notes)
    if cfg["run"]["tmin"] is not None: FR = [f for f in FR if f[0] >= cfg["run"]["tmin"] - 1e-9]
    if cfg["run"]["tmax"] is not None: FR = [f for f in FR if f[0] <= cfg["run"]["tmax"] + 1e-9]
    if len(FR) < 2:
        sys.exit("fewer than two snapshots")
    nominal_gaps([f[0] for f in FR], cfg["run"]["nominal_dt_ns"], notes)
    print(f"{len(FR)} snapshots, t = {FR[0][0]:g} .. {FR[-1][0]:g} ns")

    model = run_parameters(run, base, cfg, notes)

    # reference mesh, first and final frames, selection
    F0 = read_frame(FR[0][1], cfg["fields"])
    FL = read_frame(FR[-1][1], cfg["fields"])
    geo = Geo(F0["pts"], F0["tri"], F0["theta"], F0["phi"])
    segL = geo.segment(FL["theta"], FL["phi"], cfg["segmentation"])
    if segL is None:
        sys.exit("segmentation of the final frame found no plateau")
    sel = select(geo, F0, FL, segL, cfg, notes)
    for s in sel:
        G_geoA[s["id"]] = geo.A[np.array(s["patch_tris"])]

    # per-frame pass (parallel)
    state = dict(geo=geo, cfg=cfg, sel=sel)
    jobs = [(i, t, p) for i, (t, p, _) in enumerate(FR)]
    if a.np > 1:
        with Pool(a.np, initializer=init_worker, initargs=(state,)) as pool:
            R = pool.map(eval_frame, jobs, chunksize=1)
    else:
        init_worker(state); R = [eval_frame(j) for j in jobs]
    bad = [r for r in R if "error" in r]
    if bad:
        sys.exit(f"{bad[0]['error']} (frame t = {bad[0]['t']})")
    print(f"frames evaluated in {time.time() - T0:.0f} s")

    # histories
    H = {"t": np.array([r["t"] for r in R])}
    for j, s in enumerate(sel):
        H[s["id"]] = {k: np.array([r["points"][j][k] for r in R]) for k in R[0]["points"][j]}
        H[s["id"] + "_patch"] = {k: np.array([r["patches"][j][k] for r in R]) for k in ("theta", "theta_sd", "phi", "phi_sd", "rhoS",
                                                                                        "rhoS_sd", "rhoG", "rhoG_sd", "rhoF", "rhoF_sd")}
        H[s["id"] + "_patch"]["slipE"] = np.array([r["patches"][j]["slipE"] for r in R])
        for h in (H[s["id"]], H[s["id"] + "_patch"]):   # accumulated slip rate: backward difference between saved snapshots
            dS = np.diff(h["slipE"]); dt_ = np.diff(H["t"])
            h["slip_rate"] = np.r_[np.nan, np.maximum(dS, 0.0) / np.maximum(dt_, 1e-12)]
            h["slip_rate_sd"] = np.zeros_like(h["slip_rate"])
        for h in (H[s["id"]], H[s["id"] + "_patch"]):   # forest density as % of its own maximum over all analysed frames
            mx = np.nanmax(h["rhoF"]) if np.any(np.isfinite(h["rhoF"])) else np.nan
            h["rhoF_pct"] = 100.0 * h["rhoF"] / mx if mx and np.isfinite(mx) and mx > 0 else np.full_like(h["rhoF"], np.nan)
            if "rhoF_sd" in h:
                h["rhoF_pct_sd"] = 100.0 * h["rhoF_sd"] / mx if mx and np.isfinite(mx) and mx > 0 else np.full_like(h["rhoF"], np.nan)
    rates = reconstruct_rates(R, sel, model)

    # accumulated slip rate maps: backward difference of slip_rate_path between consecutive saved frames, written to _frames/
    tmpd = os.path.join(out, "_frames"); os.makedirs(tmpd, exist_ok=True)
    pos, nreset = [], 0
    for n in range(len(R)):
        cur = R[n].pop("slipE_all")
        if n == 0 or cur is None or prev is None:
            rt = np.zeros(len(geo.T), np.float32) if cur is not None else None
        else:
            dS = cur - prev; nreset += int((dS < -1e-12).sum())
            rt = (np.maximum(dS, 0) / max(R[n]["t"] - R[n - 1]["t"], 1e-12)).astype(np.float32)
        if rt is not None:
            np.save(os.path.join(tmpd, f"rate_{n:05d}.npy"), rt)
            if n > 0 and (rt > 0).any():
                pos.append(rt[rt > 0][:: max(1, int((rt > 0).sum() // 20000))])
        prev = cur
    notes["slip_rate"] = dict(field=cfg["fields"].get("slip_rate_path"), negative_increments=nreset,
                              available=bool(len(pos)))

    # fixed colour scales over the whole run
    sc = cfg["scales"]; L = {}
    floor = sc["vis_floor_um2"]
    th_lo = min(r["lim"]["theta"]["min"] for r in R); th_hi = max(r["lim"]["theta"]["max"] for r in R)
    if sc["theta_deg"]:
        L["theta"] = tuple(float(x) for x in sc["theta_deg"])
    else:   # convention of grain_post.py / grain_count.py: the first frame's range +- pad (twilight)
        L["theta"] = (R[0]["lim"]["theta"]["min"] - sc["theta_pad_deg"], R[0]["lim"]["theta"]["max"] + sc["theta_pad_deg"])
    clip = {"theta": dict(range=L["theta"], data_min=th_lo, data_max=th_hi, clipped=bool(th_lo < L["theta"][0] or th_hi > L["theta"][1]))}
    for k in ("rhoS", "rhoG"):
        mx = max(r["lim"][k]["max"] for r in R)
        pmin = np.nanmin([r["lim"][k]["pos_min"] for r in R])
        L[k] = tuple(float(x) for x in sc[k + "_um2"]) if sc[k + "_um2"] else (float(max(floor, pmin)), float(max(mx, 10 * floor)))
        clip[k] = dict(range=L[k], data_max=mx, data_pos_min=float(pmin),
                       below_display_floor_cells=int(sum(r["lim"][k]["n_zero"] + r["lim"][k]["n_neg"] for r in R)),
                       negative_cells_total=int(sum(r["lim"][k]["n_neg"] for r in R)),
                       clipped_above=bool(mx > L[k][1]))
    L["xlim"] = (geo.P[:, 0].min(), geo.P[:, 0].max()); L["ylim"] = (geo.P[:, 1].min(), geo.P[:, 1].max())
    if sc.get("slip_rate_per_ns"):
        L["slip_rate"] = tuple(float(x) for x in sc["slip_rate_per_ns"])
    elif pos:   # same convention as animate.py: 5th and 99.8th percentiles of the positive values over the run
        a_ = np.concatenate(pos); L["slip_rate"] = (float(max(np.percentile(a_, 5), 1e-6)), float(max(np.percentile(a_, 99.8), 1e-5)))
    else:
        L["slip_rate"] = (1e-6, 1e-5)
    clip["slip_rate"] = dict(range=L["slip_rate"], negative_increments_shown_as_zero=nreset)
    notes["colour_scales"] = clip

    # ---------------------------------------------------------------- CSVs
    import csv
    with open(os.path.join(out, "selection.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["point_id", "kind", "label", "X_um", "Y_um", "patch_radius_um", "patch_area_um2", "patch_triangles",
                    "final_boundary_dist_um", "final_grain_id(label only)", "rationale"])
        for s in sel:
            w.writerow([s["id"], s["kind"], s["label"], f"{s['x']:.6f}", f"{s['y']:.6f}", f"{s['patch_radius_um']:.4f}",
                        f"{s['patch_area_um2']:.5f}", len(s["patch_tris"]), f"{s['final_boundary_dist_um']:.4f}", s["final_grain"], s["rationale"]])
    with open(os.path.join(out, "point_history.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["t_ns", "point_id", "X_um", "Y_um", "theta_t_deg(raw stored)", "rhoS_um-2", "rhoG_um-2",
                    "rhoF_um-2(rhoS+rhoG)", "rhoF_pct_of_own_max", "slip_rate_per_ns(since previous snapshot)", "phi",
                    "slip_path_phys", "dist_to_boundary_um(instantaneous)", "grain_id(instantaneous label)", "grain_mean_theta_deg"])
        for n, r in enumerate(R):
            for j, s in enumerate(sel):
                p = r["points"][j]
                w.writerow([f"{r['t']:.6g}", s["id"], f"{s['x']:.6f}", f"{s['y']:.6f}", f"{p['theta']:.6f}", f"{p['rhoS']:.6g}",
                            f"{p['rhoG']:.6g}", f"{p['rhoF']:.6g}", f"{H[s['id']]['rhoF_pct'][n]:.4f}", f"{H[s['id']]['slip_rate'][n]:.6g}", f"{p['phi']:.6f}", f"{p['slip']:.6g}", f"{p['bdist']:.5f}", p["grain"], f"{p['grain_theta']:.4f}"])
    with open(os.path.join(out, "patch_history.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["t_ns", "patch_id", "patch_area_um2(fixed reference)", "theta_t_mean_deg", "theta_t_sd_deg", "rhoS_mean", "rhoS_sd",
                    "rhoG_mean", "rhoG_sd", "rhoF_mean(rhoS+rhoG)", "rhoF_sd", "rhoF_pct_of_own_max", "slip_rate_mean_per_ns", "phi_mean", "phi_sd", "cum_storage_est", "cum_dynrec_est", "cum_sweep_est",
                    "cum_net_pred_est", "cum_net_observed", "rates_source"])
        for n, r in enumerate(R):
            for j, s in enumerate(sel):
                p = r["patches"][j]; q = rates[j][n]
                w.writerow([f"{r['t']:.6g}", s["id"], f"{s['patch_area_um2']:.6f}", f"{p['theta']:.6f}", f"{p['theta_sd']:.4g}",
                            f"{p['rhoS']:.6g}", f"{p['rhoS_sd']:.4g}", f"{p['rhoG']:.6g}", f"{p['rhoG_sd']:.4g}", f"{p['rhoF']:.6g}",
                            f"{p['rhoF_sd']:.4g}", f"{H[s['id'] + '_patch']['rhoF_pct'][n]:.4f}", f"{H[s['id'] + '_patch']['slip_rate'][n]:.6g}", f"{p['phi']:.6f}",
                            f"{p['phi_sd']:.4g}", f"{q['store']:.6g}", f"{q['dyn']:.6g}", f"{q['sweep']:.6g}", f"{q['pred']:.6g}",
                            f"{q['obs']:.6g}", "saved rate fields" if q["exact"] else "snapshot estimate"])

    # ---------------------------------------------------------------- static figures
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    col = style_for(sel)
    tl = cfg["timeline"]

    def markers(ax):
        for m in tl["markers_ns"]:
            if H["t"][0] <= m <= H["t"][-1]:
                ax.axvline(m, color="0.4", ls="--", lw=0.8)

    fig, axs = plt.subplots(7, 1, figsize=(9, 16), sharex=True)
    for s in sel:
        h, hp = H[s["id"]], H[s["id"] + "_patch"]
        for ax, k in zip(axs[:6], ("theta", "rhoS", "rhoG", "rhoF", "rhoF_pct", "slip_rate")):
            ax.plot(H["t"], h[k], "-", color=col[s["id"]], lw=1.0, label=f"{s['id']} point")
            ax.plot(H["t"], hp[k], ":", color=col[s["id"]], lw=1.2, label=f"{s['id']} patch mean")
            ax.fill_between(H["t"], hp[k] - hp[k + "_sd"], hp[k] + hp[k + "_sd"], color=col[s["id"]], alpha=0.12, lw=0)
        axs[6].plot(H["t"], h["bdist"], "-", color=col[s["id"]], lw=1.0)
    for ax, lab in zip(axs, (r"$\theta_t$ (deg, stored)", r"$\rho_S$ ($\mu$m$^{-2}$)", r"$\rho_G$ ($\mu$m$^{-2}$)",
                             r"$\rho_F=\rho_S+\rho_G$ ($\mu$m$^{-2}$)", "$\\rho_F$ (% of its own\nmax over all t)",
                             "accumulated slip rate\n(1/ns, since the\nprevious snapshot)",
                             "distance to the\ninstantaneous boundary ($\\mu$m)")):
        ax.set_ylabel(lab, fontsize=9); markers(ax)
    axs[1].set_yscale("log"); axs[2].set_yscale("log"); axs[3].set_yscale("log"); axs[4].set_ylim(0, 105); axs[5].set_yscale("log")
    axs[0].legend(fontsize=6, ncol=4, frameon=False)
    axs[-1].set_xlabel("t (ns)")
    fig.suptitle("Tracked material points (solid) and fixed reference patches (dotted; band = area-weighted spatial s.d., "
                 "not a confidence interval)", fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(out, "histories." + ext), dpi=200)
    plt.close(fig)

    for j, s in enumerate(sel):
        q = rates[j]
        tt = np.array([r["t"] for r in q])
        fig, ax = plt.subplots(figsize=(7.5, 4.2))
        ax.plot(tt, [r["obs"] for r in q], "k-", lw=1.6, label=r"observed $\Delta\langle\rho_S\rangle_P$")
        ax.plot(tt, [r["store"] for r in q], "-", color="#2ca02c", label="storage (cumulative)")
        ax.plot(tt, [-r["dyn"] for r in q], "-", color="#1f77b4", label="- dynamic recovery")
        ax.plot(tt, [-r["sweep"] for r in q], "-", color="#d62728", label="- sweeping recovery")
        ax.plot(tt, [r["pred"] for r in q], "--", color="0.3", label="predicted net (sum)")
        markers(ax)
        src = "saved rate fields" if any(r["exact"] for r in q) else "estimated from snapshots (see summary.txt)"
        ax.set_title(f"{s['id']} ({s['label']}): patch SSD budget, {src}", fontsize=9)
        ax.set_xlabel("t (ns)"); ax.set_ylabel(r"$\mu$m$^{-2}$"); ax.legend(fontsize=7, frameon=False)
        fig.tight_layout()
        for ext in ("pdf", "png"):
            fig.savefig(os.path.join(out, f"recovery_{s['id']}." + ext), dpi=200)
        plt.close(fig)

    # snapshot grid: requested times (nearest saved frame within half a spacing) + first frame of each new grain's plateau
    req = []
    for tq in tl["snapshot_times_ns"]:
        k = int(np.argmin(np.abs(H["t"] - tq)))
        if abs(H["t"][k] - tq) <= 0.5 * cfg["run"]["nominal_dt_ns"] + 1e-9:
            req.append(k)
        else:
            notes["warnings"].append(f"snapshot time {tq} ns requested but not available (nearest {H['t'][k]:g} ns)")
    ev = []
    for s in sel:
        if s["kind"] != "new-interior":
            continue
        th = H[s["id"]]["theta"]
        k = np.flatnonzero(np.abs(th - th[0]) >= cfg["selection"]["new_grain_dtheta_deg"])
        if len(k):
            ev.append(int(k[0])); notes.setdefault("formation_events", []).append(
                f"{s['id']}: theta_t first differs from its first-frame value by >= {cfg['selection']['new_grain_dtheta_deg']} deg at t = {H['t'][k[0]]:g} ns")
    idx = sorted(set(req + ev))
    notes["snapshot_frames"] = [float(H["t"][k]) for k in idx]
    G.update(dict(geo=geo, cfg=cfg, sel=sel, R=R, H={"t": H["t"], **{s["id"]: H[s["id"]] for s in sel}}, L=L))
    from scipy.spatial import cKDTree
    G["kd"] = cKDTree(geo.P); G["tmpd"] = tmpd
    snaps = [render_frame((k, os.path.join(tmpd, f"snap_{k:05d}.png"))) for k in idx]
    from PIL import Image
    ims = [Image.open(p).convert("RGB") for p in snaps]
    W_, H_ = ims[0].size
    sheet = Image.new("RGB", (W_, H_ * len(ims)), "white")
    for i, im in enumerate(ims):
        sheet.paste(im, (0, i * H_))
    sheet.save(os.path.join(out, "snapshots.png"))
    ims[0].save(os.path.join(out, "snapshots.pdf"), save_all=True, append_images=ims[1:], resolution=cfg["render"]["dpi"])

    # relaxation diagnostics
    lg = os.path.join(run, base + "_log.csv")
    if os.path.exists(lg):
        raw = np.genfromtxt(lg, delimiter=",", names=True)
        d = {k: np.atleast_1d(raw[k]).astype(float) for k in raw.dtype.names}
        tt = d["t_ns"]; KE = d.get("KE_trans", np.zeros_like(tt)) + d.get("KE_rot", np.zeros_like(tt))
        fig, axs = plt.subplots(3, 1, figsize=(8, 8), sharex=True)
        axs[0].semilogy(tt, np.maximum(KE, 1e-12), "k-", lw=1)
        axs[0].set_ylabel("kinetic energy (GPa um^2)")
        if "W_plate" in d and "Q_heat" in d:
            W = np.where(d["W_plate"] != 0, d["W_plate"], np.nan)
            axs[1].plot(tt, d["Q_heat"] / W, "r-", lw=1); axs[1].set_ylabel("Q_heat / W_plate")
            for tq in (13.0, 30.0):
                k = int(np.argmin(np.abs(tt - tq)))
                if abs(tt[k] - tq) < 0.06:
                    notes[f"energy_at_{tq:g}ns"] = dict(KE=float(KE[k]), Q_over_W=float(d["Q_heat"][k] / W[k]) if np.isfinite(W[k]) else None)
        else:
            notes["warnings"].append("the log has no W_plate/Q_heat columns: no heat-conversion check")
        st = [(r["t"], r["stress"].get("sig_rms", np.nan), r["stress"].get("sig_p99", np.nan)) for r in R]
        axs[2].plot([a_[0] for a_ in st], [a_[1] for a_ in st], "b-", lw=1, label="area RMS of |sigma|_F")
        axs[2].plot([a_[0] for a_ in st], [a_[2] for a_ in st], "b:", lw=1, label="99th percentile of |sigma|_F")
        axs[2].set_ylabel("stress (GPa)"); axs[2].legend(fontsize=8, frameon=False); axs[2].set_xlabel("t (ns)")
        for ax in axs:
            markers(ax)
        fig.suptitle("Relaxation diagnostics (energy ledger <base>_log.csv; stress from the saved snapshots)", fontsize=9)
        fig.tight_layout(rect=(0, 0, 1, 0.97))
        for ext in ("pdf", "png"):
            fig.savefig(os.path.join(out, "relaxation." + ext), dpi=200)
        plt.close(fig)
        notes["energy_ledger"] = lg
    else:
        notes["warnings"].append("no <base>_log.csv: no energy ledger, relaxation documented by the stress diagnostic only")

    # ---------------------------------------------------------------- animation
    if cfg["render"]["mp4"]:
        variants = [False] + ([True] if cfg["render"]["deformed"] else [])
        for dfm in variants:
            cfg["render"]["_deformed"] = dfm
            G["cfg"] = cfg
            pngs = [(n, os.path.join(tmpd, f"{'d' if dfm else 'r'}_{n:05d}.png")) for n in range(len(R))]
            if a.np > 1:
                with Pool(a.np, initializer=init_worker, initargs=(dict(G),)) as pool:
                    pool.map(render_frame, pngs, chunksize=1)
            else:
                for p in pngs:
                    render_frame(p)
            mp4 = os.path.join(out, "grain_recovery_deformed.mp4" if dfm else "grain_recovery.mp4")
            if shutil.which("ffmpeg"):
                cmd = ["ffmpeg", "-y", "-loglevel", "error", "-framerate", str(cfg["render"]["fps"]), "-i",
                       os.path.join(tmpd, f"{'d' if dfm else 'r'}_%05d.png"), "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2",
                       "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "20", mp4]
                subprocess.run(cmd, check=True)
            else:
                notes["warnings"].append("ffmpeg not found: frames left in _frames/, no MP4")
        cfg["render"].pop("_deformed", None)
        if shutil.which("ffmpeg"):
            for f in os.listdir(tmpd):
                if not f.startswith("snap_"):
                    os.remove(os.path.join(tmpd, f))

    # ---------------------------------------------------------------- summary
    write_summary(out, cfg, notes, FR, R, H, sel, rates, model)
    rc = json.loads(json.dumps(cfg))
    rc["selection"]["mode"] = "manual"
    rc["selection"]["points"] = [dict(x=s["x"], y=s["y"], kind=s["kind"], label=s["label"], rationale=s["rationale"],
                                      patch_radius_um=s["patch_radius_um"]) for s in sel]
    json.dump(rc, open(os.path.join(out, "resolved_config.json"), "w"), indent=1)
    print(f"done in {time.time() - T0:.0f} s -> {out}")


def write_summary(out, cfg, notes, FR, R, H, sel, rates, model):
    L = []
    w = L.append
    w("Grain formation, recovery and retention: diagnostic summary")
    w("=" * 60)
    w(f"run: {cfg['run']['dir']}  base: {cfg['run']['base']}  analysed: {notes['started']}")
    w(f"code version recorded by the run: {notes.get('code_version')}")
    w("")
    w("SNAPSHOTS")
    w(f"  analysed: {len(FR)} frames, t = {FR[0][0]:g} .. {FR[-1][0]:g} ns ({sum(1 for f in FR if f[2] == 'relaxed')} relaxed frame at t = 0)")
    w(f"  .pvd entries: {notes.get('pvd_entries')}; listed but missing on disk: {len(notes['listed_but_missing'])}"
      + (f" ({', '.join(f'{t:g}' for t, _ in notes['listed_but_missing'][:20])}{' ...' if len(notes['listed_but_missing']) > 20 else ''} ns)" if notes['listed_but_missing'] else ""))
    w(f"  duplicate times dropped: {len(notes['duplicates_dropped'])}")
    ab = notes.get("nominal_times_absent", [])
    w(f"  nominal {cfg['run']['nominal_dt_ns']} ns grid points without a frame: {len(ab)}" + (f" (first: {', '.join(f'{x:g}' for x in ab[:15])})" if ab else ""))
    if notes.get("off_grid_times"):
        w(f"  frames off the nominal grid: {', '.join(f'{x:g}' for x in notes['off_grid_times'])} ns")
    w(f"  expected 301 frames for 0..30 ns at 0.1 ns only if both end points are present; found {len(FR)}")
    w("")
    w("MATERIAL CORRESPONDENCE")
    w("  every frame has the same points and triangles as the first one (checked by hash): the mesh is a fixed material")
    w("  (reference) mesh, so node and element indices identify the same material in every frame, also for MPI runs (rank 0")
    w("  writes the global mesh). Points: fixed (triangle, barycentric) positions; P1 fields (theta_t, phi) interpolated")
    w("  linearly, P0 fields (rhoS, rhoG) taken from the containing triangle. Patches: fixed triangle sets, weights = reference")
    w("  areas (planar model); P1 fields enter a patch mean through their exact triangle means.")
    w("")
    s_ = cfg["segmentation"]
    w("SEGMENTATION (grain_count.segment, unchanged)")
    w(f"  field theta_lattice_deg (stored total orientation, no wrapping or symmetry reduction applied to it);")
    w(f"  plateau gradient threshold {s_['grad_deg_per_um']} deg/um, minimum plateau area {s_['minarea_um2']} um^2, merge {s_['merge_deg']} deg;")
    w(f"  symmetry period {s_['sym_deg']} deg used only for the misorientation reported by grain_count (not for display or tracking).")
    w(f"  GND density is NOT used for the segmentation. Final frame: {notes['final_grains']} grains, {notes['new_grain_candidates']} new-grain candidates")
    w(f"  (criterion: plateau rotated >= {cfg['selection']['new_grain_dtheta_deg']} deg since the first frame, or >= that from the remnant of its initial grain).")
    w("")
    w("SELECTED POINTS AND PATCHES (fixed in reference coordinates; final grain ids are labels, not tracked lineages)")
    for s in sel:
        w(f"  {s['id']} [{s['kind']}] X = ({s['x']:.3f}, {s['y']:.3f}) um, patch r = {s['patch_radius_um']:.3f} um "
          f"({len(s['patch_tris'])} triangles, {s['patch_area_um2']:.4f} um^2), final boundary/edge distance {s['final_boundary_dist_um']:.3f} um")
        w(f"      {s['rationale']}")
    w("")
    w("DISPLAY")
    for k, v in notes["colour_scales"].items():
        w(f"  {k}: fixed range {tuple(round(x, 6) for x in v['range'])} over the whole run; " + ", ".join(f"{a}={b}" for a, b in v.items() if a != "range"))
    w(f"  density display floor {cfg['scales']['vis_floor_um2']} um^-2 (display only; the solver floors are rhoS >= {model['rho_floor']} and")
    w("  rho_m >= 1e-3, unrelated); CSVs carry the raw values. Negative densities are counted above, never made positive.")
    w("  Playback advances one saved frame per video frame; no interpolated frames.")
    sr = notes.get("slip_rate", {})
    w(f"  accumulated slip rate: backward difference of the saved accumulated slip '{sr.get('field')}' (sum over systems of |dgamma|,")
    w("  as in the 'slip_acc' animation) between consecutive snapshots, divided by their time difference and shown at the later one:")
    w("  an average over the snapshot interval, not the instantaneous rate; no value at the first frame;")
    w(f"  {sr.get('negative_increments', 0)} element increments were negative (a reset of the slip path) and are shown as 0.")
    w("")
    w("SSD LAW AND RATE RECONSTRUCTION")
    w(f"  parameters: K = {model['K']}, d = {model['d_um']} um, b = {model['b_um']} um, C_D = {model['C_D']}, C_A = {model['C_A_um']} um, "
      f"recovery mode -recov {model['recov']}, floor {model['rho_floor']} um^-2")
    for k, v in notes["model_parameter_sources"].items():
        w(f"    {k}: {v}")
    w("  law (per solver step): drho = K/b sqrt(rho) dGamma - 2d/b rho dGamma - C_D <dphi>+ tanh(C_A^2 |grad theta_t|^2) rho (mode 2),")
    w("  storage capped so that Z drho_store <= sum |tau_a dgamma_a| (not reconstructable from snapshots), floor applied after.")
    exact = any(r["exact"] for q in rates for r in q)
    if exact:
        w("  saved rate fields were used (rectangle rule over the snapshot interval).")
    else:
        w("  no saved rate fields: storage and dynamic recovery are ESTIMATES from the exact cumulative slip path slip_path_phys")
        w("  (Gamma) with the density averaged over each snapshot interval; sweeping recovery is an ESTIMATE from the snapshot")
        w("  change of the element-mean phi (positive part taken after differencing: episodes of increasing phi that are reversed")
        w("  within an interval are missed, so the sweep estimate is a lower bound of the per-step sum).")
    w("  Budget per patch (cumulative over the analysed interval; um^-2):")
    for j, s in enumerate(sel):
        r = rates[j][-1]
        mis = r["obs"] - r["pred"]
        w(f"    {s['id']}: observed {r['obs']:+.4g}; storage {r['store']:+.4g}, dyn. recovery {-r['dyn']:+.4g}, sweep {-r['sweep']:+.4g}; "
          f"predicted net {r['pred']:+.4g}; mismatch {mis:+.4g}")
    w("  A mismatch comes from snapshot sampling (density not constant over an interval, unresolved phi episodes), the storage")
    w("  cap and the floor; it is not an exact discrete budget.")
    w("")
    w("RESETS")
    if FR[0][2] == "relaxed":
        r0 = R[0]
        w("  t = 0 is the relaxed frame. At the start of loading the solver resets rho_S = rho_M = rho0 everywhere; a change between")
        w("  the relaxed frame and the first loaded frame that matches that reset is not physical recovery.")
    else:
        w("  no relaxed frame in the analysed set; the first frame is the first loaded snapshot.")
    w("")
    w("OBSERVATIONS (computed, not interpreted)")
    for j, s in enumerate(sel):
        h = H[s["id"]]; hp = H[s["id"] + "_patch"]
        k1 = int(np.argmax(hp["rhoS"]))
        w(f"  {s['id']} [{s['kind']}]: patch rhoS {hp['rhoS'][0]:.4g} -> max {hp['rhoS'][k1]:.4g} at {H['t'][k1]:g} ns -> {hp['rhoS'][-1]:.4g} at "
          f"{H['t'][-1]:g} ns; rhoG {hp['rhoG'][0]:.4g} -> {hp['rhoG'][-1]:.4g}; phi min {h['phi'].min():.3f}; "
          f"min distance to the instantaneous boundary {np.nanmin(h['bdist']):.3f} um")
        kf = int(np.nanargmax(hp["rhoF"]))
        w(f"      patch rhoF = rhoS + rhoG: max {hp['rhoF'][kf]:.4g} at {H['t'][kf]:g} ns; at {H['t'][-1]:g} ns {hp['rhoF'][-1]:.4g} "
          f"= {hp['rhoF_pct'][-1]:.1f}% of that maximum")
        dep = hp["rhoS"][-1] < hp["rhoS"][k1] * 0.9
        w(f"      net SSD depletion after the maximum: {'yes' if dep else 'no'} ({(1 - hp['rhoS'][-1] / hp['rhoS'][k1]) * 100:.1f}% below the maximum)")
    for e in notes.get("formation_events", []):
        w("  " + e)
    w(f"  snapshot figure frames: {', '.join(f'{x:g}' for x in notes['snapshot_frames'])} ns")
    for tq in (13.0, 30.0):
        if f"energy_at_{tq:g}ns" in notes:
            e = notes[f"energy_at_{tq:g}ns"]
            w(f"  at {tq:g} ns: kinetic energy {e['KE']:.4g} GPa um^2, Q_heat/W_plate = {e['Q_over_W']}")
    w("")
    w("CAVEATS")
    w("  Events faster than the snapshot spacing are unresolved. Retention to the last frame shows persistence over the simulated")
    w("  interval only. Agreement of GND bands with orientation boundaries is internal consistency of coupled model fields, not")
    w("  experimental validation. The saved rhoG is the measure of the code version that wrote it (versions before 7.5:")
    w("  |F_p c - G_R|/b; 7.5: Cermelli-Gurtin |c/J_p|/b without reference).")
    if notes["warnings"]:
        w(""); w("WARNINGS")
        for x in notes["warnings"]:
            w("  " + x)
    open(os.path.join(out, "summary.txt"), "w").write("\n".join(L) + "\n")


if __name__ == "__main__":
    main()
