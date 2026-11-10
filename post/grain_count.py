#!/usr/bin/env python3
"""grain_count.py  --  count distinct grains and their neighbour misorientations in a solver .vtu frame.

usage:  python3 grain_count.py frame.vtu [-o prefix] [--grad 2] [--merge 2] [--minarea 0.1] [--sym 90] [--scan]

Grains are PLATEAUS of the lattice orientation (point field theta_lattice_deg).  Smooth ramps between plateaus
(e.g. 8 deg over ~1 um) are not grains: they are the boundary between the two plateaus they connect.

Method (needs only numpy, scipy, matplotlib; the VTU is read directly, no meshio/vtk):
  1. element gradient |grad theta| (deg/um) on the triangles (P1 field);
  2. plateau seeds = connected sets of triangles with |grad theta| < --grad (deg/um) and area >= --minarea (um^2);
  3. every other triangle (ramp, sharp core, small sliver) goes to the nearest seed along the mesh (graph distance);
  4. neighbouring plateaus whose mean orientations differ by < --merge (deg) are merged (agglomerative, closest first;
     a merged grain keeps every member within --merge of its mean), so a smoothly bent grain is not split into pieces;
  5. grains = the merged plateaus; neighbours = grains whose territories share mesh edges;
     misorientation = difference of the plateau-mean orientations; width of the boundary = median over the shared edges
     of (distance to the nearest plateau of grain i + distance to the nearest plateau of grain j), i.e. the width of the
     transition between the two plateaus (a sharp KWC core gives ~ one mesh size, your 8 deg / 1 um ramp gives ~ 1).
Misorientation is reported twice: raw (theta_i - theta_j, as the solver sees it) and reduced by the crystal symmetry
period --sym (90 deg = cubic, rotation about [001]; 360 = no symmetry).
--probe x,y,r explains what the segmentation did at one place (use it for a grain that seems to be missed).
--scan prints the grain count for a grid of --grad / --merge values: use it to see how robust the count is.
Also prints the grain-boundary LENGTH fraction per misorientation bin ([2,5), [5,15), >=15 deg; --lo/--mid/--hi) -> <prefix>_misorient_hist.csv.
Outputs: <prefix>_grains.csv, <prefix>_boundaries.csv, <prefix>_misorient_hist.csv and three figures <prefix>_theta.png, <prefix>_territory.png, <prefix>_misorientation.png.
"""
import argparse, base64, re, sys
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components, dijkstra
from scipy.spatial import cKDTree

DT = {"Float32": "<f4", "Float64": "<f8", "Int32": "<i4", "Int64": "<i8", "UInt8": "u1", "UInt32": "<u4", "UInt64": "<u8", "Int8": "i1"}


def read_vtu(fn):
    """Minimal reader for the solver's inline-base64 .vtu (uncompressed, UInt32 block header)."""
    txt = open(fn, "r", errors="ignore").read()
    out = {}
    for m in re.finditer(r"<DataArray([^>]*)>(.*?)</DataArray>", txt, re.S):
        at, body = m.group(1), m.group(2).strip()
        name = re.search(r'Name="([^"]*)"', at)
        typ = re.search(r'type="([^"]*)"', at).group(1)
        nc = re.search(r'NumberOfComponents="(\d+)"', at)
        if 'format="binary"' not in at:
            continue
        raw = base64.b64decode("".join(body.split()))
        n = int(np.frombuffer(raw[:4], "<u4")[0])
        a = np.frombuffer(raw[4:4 + n], DT[typ])
        if nc and int(nc.group(1)) > 1:
            a = a.reshape(-1, int(nc.group(1)))
        out[name.group(1) if name else "points"] = a
    return out


def load(fn):
    d = read_vtu(fn)
    pts = d["points"][:, :2].astype(float)
    conn, off, typ = d["connectivity"], d["offsets"], d["types"]
    start = np.r_[0, off[:-1]]
    tri = np.array([conn[s:s + 3] for s, t in zip(start, typ) if t == 5])
    istri = typ == 5
    label = d["Label"][istri] if "Label" in d else np.zeros(len(tri), int)
    th = d["theta_lattice_deg"].astype(float).ravel()
    phi = d["phi"].astype(float).ravel() if "phi" in d else np.zeros(len(pts))
    return pts, tri, th, phi, label


class Mesh:
    def __init__(self, pts, tri, th, phi):
        self.P, self.T, self.th, self.phi = pts, tri, th, phi
        p0, p1, p2 = pts[tri[:, 0]], pts[tri[:, 1]], pts[tri[:, 2]]
        e1, e2 = p1 - p0, p2 - p0
        det = e1[:, 0] * e2[:, 1] - e1[:, 1] * e2[:, 0]
        d1, d2 = th[tri[:, 1]] - th[tri[:, 0]], th[tri[:, 2]] - th[tri[:, 0]]
        self.g = np.hypot((d1 * e2[:, 1] - d2 * e1[:, 1]) / det, (-d1 * e2[:, 0] + d2 * e1[:, 0]) / det)  # deg/um
        self.A = np.abs(det) / 2
        self.cen = (p0 + p1 + p2) / 3
        self.tm = th[tri].mean(1)
        self.h = np.sqrt(2 * self.A.mean())
        # interior edges: pairs of triangles sharing an edge
        nt = len(tri)
        ed = np.concatenate([tri[:, [0, 1]], tri[:, [1, 2]], tri[:, [2, 0]]])
        tid = np.tile(np.arange(nt), 3)
        ed.sort(axis=1)
        key = ed[:, 0].astype(np.int64) * len(pts) + ed[:, 1]
        o = np.argsort(key, kind="stable")
        k, ed, tid = key[o], ed[o], tid[o]
        same = np.flatnonzero(k[1:] == k[:-1])
        self.pa, self.pb = tid[same], tid[same + 1]          # triangle pair
        self.en = ed[same]                                    # edge nodes
        self.elen = np.hypot(*(pts[self.en[:, 0]] - pts[self.en[:, 1]]).T)
        self.emid = 0.5 * (pts[self.en[:, 0]] + pts[self.en[:, 1]])
        w = np.hypot(*(self.cen[self.pa] - self.cen[self.pb]).T)
        self.graph = coo_matrix((np.r_[w, w], (np.r_[self.pa, self.pb], np.r_[self.pb, self.pa])), shape=(nt, nt)).tocsr()


def segment(M, grad, minarea, merge):
    nt = len(M.T)
    ok = M.g < grad
    sel = ok[M.pa] & ok[M.pb]
    n, lab = connected_components(coo_matrix((np.ones(sel.sum()), (M.pa[sel], M.pb[sel])), shape=(nt, nt)), directed=False)
    ar = np.bincount(lab, weights=M.A * ok, minlength=n)
    keep = np.flatnonzero(ar >= minarea)
    if len(keep) == 0:
        return None
    sub = -np.ones(nt, int)                                   # plateau id per triangle (-1 = not plateau)
    for s, c in enumerate(keep):
        sub[(lab == c) & ok] = s
    S = len(keep)
    # territories: nearest plateau triangle along the mesh
    src = np.flatnonzero(sub >= 0)
    dist, pred, srcs = dijkstra(M.graph, directed=False, indices=src, min_only=True, return_predecessors=True)
    terr = np.where(srcs >= 0, sub[np.maximum(srcs, 0)], -1)
    # plateau means
    pa_ = np.bincount(sub[src], weights=M.A[src], minlength=S)
    mu = np.bincount(sub[src], weights=(M.A * M.tm)[src], minlength=S) / pa_
    # adjacency of territories
    ta, tb = terr[M.pa], terr[M.pb]
    dif = (ta != tb) & (ta >= 0) & (tb >= 0)
    adj = {tuple(sorted(p)) for p in zip(ta[dif], tb[dif])}
    # agglomerative merge of plateaus that differ by < merge
    grp = {s: dict(mem=[s], area=pa_[s], mu=mu[s]) for s in range(S)}
    nb = {s: set() for s in range(S)}
    for a, b in adj:
        nb[a].add(b); nb[b].add(a)
    while True:
        best = None
        for a in nb:
            for b in nb[a]:
                if b <= a:
                    continue
                d = abs(grp[a]["mu"] - grp[b]["mu"])
                if d >= merge or (best is not None and d >= best[0]):
                    continue
                mem = grp[a]["mem"] + grp[b]["mem"]
                ar2 = np.array([pa_[m] for m in mem]); mu2 = np.array([mu[m] for m in mem])
                m2 = (ar2 * mu2).sum() / ar2.sum()
                if np.abs(mu2 - m2).max() <= merge:
                    best = (d, a, b, mem, ar2.sum(), m2)
        if best is None:
            break
        _, a, b, mem, area, m2 = best
        grp[a] = dict(mem=mem, area=area, mu=m2)
        nb[a] |= nb[b]; nb[a].discard(a); nb[a].discard(b)
        for c in nb[b]:
            nb[c].discard(b)
            if c != a:
                nb[c].add(a)
        del grp[b], nb[b]
    gid = {g: i for i, g in enumerate(sorted(grp, key=lambda g: -grp[g]["area"]))}
    remap = np.empty(S, int)
    for g, d in grp.items():
        for m in d["mem"]:
            remap[m] = gid[g]
    G = len(grp)
    gt = np.where(terr >= 0, remap[np.maximum(terr, 0)], -1)             # territory (grain id) of every triangle
    gp = np.where(sub >= 0, remap[np.maximum(sub, 0)], -1)               # plateau triangles of the final grains
    mug = np.zeros(G);
    for g, d in grp.items():
        mug[gid[g]] = d["mu"]
    return dict(G=G, terr=gt, plat=gp, mu=mug)


def analyse(M, seg, label, sym):
    G, terr, plat, mu = seg["G"], seg["terr"], seg["plat"], seg["mu"]
    rows = []
    for g in range(G):
        p = plat == g
        t = terr == g
        a = M.A[p]
        c = (M.cen[p] * a[:, None]).sum(0) / a.sum()
        sd = np.sqrt((a * (M.tm[p] - mu[g]) ** 2).sum() / a.sum())
        lb = np.bincount(label[p], weights=a)
        rows.append(dict(id=g, plateau_area=a.sum(), territory_area=M.A[t].sum(), cx=c[0], cy=c[1], theta_deg=mu[g],
                         theta_std_deg=sd, theta_sym_deg=((mu[g] + sym / 2) % sym) - sym / 2 if sym < 360 else mu[g],
                         orig_label=int(np.argmax(lb)) if len(lb) else -1))
    # boundaries
    ta, tb = terr[M.pa], terr[M.pb]
    dif = ta != tb
    pairs = {}
    for e in np.flatnonzero(dif):
        pairs.setdefault((min(ta[e], tb[e]), max(ta[e], tb[e])), []).append(e)
    trees = [cKDTree(M.cen[plat == g]) for g in range(G)]
    bnd = []
    for (i, j), es in sorted(pairs.items()):
        es = np.array(es)
        L = M.elen[es].sum()
        di = trees[i].query(M.emid[es])[0]; dj = trees[j].query(M.emid[es])[0]
        width = np.median(di + dj)
        d = mu[j] - mu[i]
        ds = ((d + sym / 2) % sym) - sym / 2 if sym < 360 else d
        bnd.append(dict(i=i, j=j, misorient_deg=abs(d), misorient_sym_deg=abs(ds), length_um=L, width_um=width,
                        phi_at_edge=M.phi[M.en[es].ravel()].mean(), edges=es))
    return rows, bnd


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("vtu"); ap.add_argument("-o", "--out", default=None)
    ap.add_argument("--grad", type=float, default=2.0, help="plateau gradient threshold, deg/um (2)")
    ap.add_argument("--merge", type=float, default=2.0, help="merge plateaus closer than this, deg (2)")
    ap.add_argument("--minarea", type=float, default=0.1, help="smallest plateau, um^2 (0.1)")
    ap.add_argument("--sym", type=float, default=90.0, help="crystal period for the reduced misorientation, deg (90 cubic, 360 none)")
    ap.add_argument("--lo", type=float, default=2.0, help="lower bin edge, deg (2)")
    ap.add_argument("--mid", type=float, default=5.0, help="middle bin edge, deg (5; bins [lo,mid) [mid,hi) >=hi)")
    ap.add_argument("--hi", type=float, default=15.0, help="upper bin edge, deg (15)")
    ap.add_argument("--scan", action="store_true", help="grain count vs --grad and --merge")
    ap.add_argument("--probe", default=None, help="x,y,r (um): why is the grain at (x,y) found / missed? prints theta, |grad theta| and the plateau/territory ids within r")
    a = ap.parse_args()
    report(a, a.out or re.sub(r"\.vtu$", "", a.vtu) + "_grains")


def bins_of(a):
    return [("<%g" % a.lo, 0, a.lo), (f"[{a.lo:g},{a.mid:g})", a.lo, a.mid), (f"[{a.mid:g},{a.hi:g})", a.mid, a.hi), (f">={a.hi:g}", a.hi, 1e9)]


def analyse_frame(fn, a):
    """Load + segment + analyse one frame: (M, seg, rows, bnd) or None if no plateau was found."""
    pts, tri, th, phi, label = load(fn)
    M = Mesh(pts, tri, th, phi)
    seg = segment(M, a.grad, a.minarea, a.merge)
    if seg is None:
        return None
    rows, bnd = analyse(M, seg, label, a.sym)
    return M, seg, rows, bnd


def summary_from(M, seg, rows, bnd, a):
    tot = sum(b["length_um"] for b in bnd) or 1e-30
    out = dict(n_grains=seg["G"], n_ge1um2=sum(r["plateau_area"] >= 1.0 for r in rows), n_boundaries=len(bnd),
               boundary_length_um=tot if bnd else 0.0, plateau_area_frac=float(M.A[seg["plat"] >= 0].sum() / M.A.sum()),
               mean_grain_area_um2=float(np.mean([r["plateau_area"] for r in rows])))
    for tag, key in (("raw", "misorient_deg"), ("sym", "misorient_sym_deg")):
        for name, lo, hi in bins_of(a):
            out[f"pct_{tag}_{name}"] = 100 * sum(b["length_um"] for b in bnd if lo <= b[key] < hi) / tot
        out[f"mean_{tag}_deg"] = sum(b[key] * b["length_um"] for b in bnd) / tot
    return out


def summary(fn, a):
    """One row of numbers for a frame."""
    r = analyse_frame(fn, a)
    return None if r is None else summary_from(*r, a)


def report(a, pfx):
    pts, tri, th, phi, label = load(a.vtu)
    M = Mesh(pts, tri, th, phi)
    print(f"{a.vtu}: {len(pts)} nodes, {len(tri)} triangles, domain {np.ptp(pts[:,0]):.2f} x {np.ptp(pts[:,1]):.2f} um, mesh size {M.h:.3f} um, "
          f"theta_lattice {th.min():.1f} .. {th.max():.1f} deg")
    if a.scan:
        gs, ms = [1, 2, 3, 5, 8], [1, 2, 3, 5]
        print("\ngrain count (rows: --grad deg/um, columns: --merge deg), --minarea %.3g" % a.minarea)
        print("grad\\merge " + "".join(f"{m:>6g}" for m in ms))
        for g in gs:
            r = []
            for m in ms:
                s = segment(M, g, a.minarea, m); r.append(0 if s is None else s["G"])
            print(f"{g:>10g} " + "".join(f"{x:>6d}" for x in r))
    seg = segment(M, a.grad, a.minarea, a.merge)
    if seg is None:
        sys.exit("no plateau found: increase --grad or decrease --minarea")
    rows, bnd = analyse(M, seg, label, a.sym)
    G = seg["G"]
    print(f"\n{G} grains  (--grad {a.grad:g} deg/um, --merge {a.merge:g} deg, --minarea {a.minarea:g} um^2),  {len(bnd)} neighbour pairs")
    big = sum(r["plateau_area"] >= 1.0 for r in rows)
    print(f"of which {big} have a plateau area >= 1 um^2 (the rest are small plateaus, usually patches of ramps and corners)")
    print(f"plateaus cover {100*M.A[seg['plat']>=0].sum()/M.A.sum():.1f} % of the area, ramps/cores {100*M.A[seg['plat']<0].sum()/M.A.sum():.1f} %\n")
    print(" id  area[um2]  centre(x,y)[um]   theta[deg]  std   theta mod %g  orig  neighbours: id:misorient(raw/sym, deg) width[um]" % a.sym)
    nbs = {g: [] for g in range(G)}
    for b in bnd:
        nbs[b["i"]].append((b["j"], b)); nbs[b["j"]].append((b["i"], b))
    for r in rows:
        s = "  ".join(f"{j}:{b['misorient_deg']:.1f}/{b['misorient_sym_deg']:.1f} w{b['width_um']:.2f}" for j, b in sorted(nbs[r["id"]]))
        print(f"{r['id']:>3d} {r['plateau_area']:>9.2f}  ({r['cx']:5.2f},{r['cy']:5.2f})  {r['theta_deg']:>9.2f} {r['theta_std_deg']:>5.2f} {r['theta_sym_deg']:>9.2f}  {r['orig_label']:>4d}  {s}")
    bins = bins_of(a)
    tot = sum(b["length_um"] for b in bnd)
    print(f"\ngrain-boundary length by misorientation (total {tot:.2f} um over {len(bnd)} boundaries):")
    print(f"{'bin [deg]':>10} | {'raw length':>10} {'%':>6} {'n':>3} | {'mod ' + format(a.sym, 'g') + ' length':>13} {'%':>6} {'n':>3}")
    hist = []
    for name, lo, hi in bins:
        r = [b["length_um"] for b in bnd if lo <= b["misorient_deg"] < hi]
        q = [b["length_um"] for b in bnd if lo <= b["misorient_sym_deg"] < hi]
        hist.append((name, sum(r), 100 * sum(r) / tot, len(r), sum(q), 100 * sum(q) / tot, len(q)))
        print(f"{name:>10} | {sum(r):>10.2f} {100*sum(r)/tot:>6.1f} {len(r):>3d} | {sum(q):>13.2f} {100*sum(q)/tot:>6.1f} {len(q):>3d}")
    with open(pfx + "_misorient_hist.csv", "w") as f:
        f.write("bin_deg,raw_length_um,raw_percent,raw_n,sym_length_um,sym_percent,sym_n\n")
        for h in hist:
            f.write(f"{h[0]},{h[1]:.4f},{h[2]:.2f},{h[3]},{h[4]:.4f},{h[5]:.2f},{h[6]}\n")
    with open(pfx + "_grains.csv", "w") as f:
        f.write("id,plateau_area_um2,territory_area_um2,cx_um,cy_um,theta_deg,theta_std_deg,theta_sym_deg,orig_label,n_neighbours\n")
        for r in rows:
            f.write(f"{r['id']},{r['plateau_area']:.4f},{r['territory_area']:.4f},{r['cx']:.3f},{r['cy']:.3f},{r['theta_deg']:.3f},"
                    f"{r['theta_std_deg']:.3f},{r['theta_sym_deg']:.3f},{r['orig_label']},{len(nbs[r['id']])}\n")
    with open(pfx + "_boundaries.csv", "w") as f:
        f.write("grain_i,grain_j,misorient_deg,misorient_sym_deg,length_um,width_um,mean_phi_at_edge\n")
        for b in bnd:
            f.write(f"{b['i']},{b['j']},{b['misorient_deg']:.3f},{b['misorient_sym_deg']:.3f},{b['length_um']:.3f},{b['width_um']:.3f},{b['phi_at_edge']:.3f}\n")
    if getattr(a, "probe", None):
        px, py, pr = [float(v) for v in a.probe.split(",")]
        m = np.hypot(M.cen[:, 0] - px, M.cen[:, 1] - py) < pr
        if m.sum() == 0:
            print("\nprobe: no triangles within r")
        else:
            A_ = M.A[m]
            print(f"\nprobe ({px:g},{py:g}) r={pr:g} um: {m.sum()} triangles, area {A_.sum():.2f} um^2")
            print(f"  theta: median {np.median(M.tm[m]):.2f}, 5-95 %: {np.percentile(M.tm[m],5):.2f} .. {np.percentile(M.tm[m],95):.2f} deg")
            print("  |grad theta| [deg/um]: 10/50/90 %% = %.2f / %.2f / %.2f" % tuple(np.percentile(M.g[m], [10, 50, 90])))
            print("  area fraction below gradient:  " + "  ".join(f"{t:g}: {100*A_[M.g[m] < t].sum()/A_.sum():.0f}%" for t in (1, 2, 3, 5, 8, 12)))
            ids, ct = np.unique(seg["plat"][m], return_counts=True)
            print("  plateau grain ids here (-1 = ramp): " + ", ".join(f"{i}: {100*c/m.sum():.0f}%" for i, c in zip(ids, ct)))
            ids, ct = np.unique(seg["terr"][m], return_counts=True)
            print("  territory grain ids here: " + ", ".join(f"{i}: {100*c/m.sum():.0f}%" for i, c in zip(ids, ct)))
    plot(M, seg, rows, bnd, pfx + ".png", a)
    print(f"\nwrote {pfx}_grains.csv, {pfx}_boundaries.csv, {pfx}_misorient_hist.csv, {pfx}_theta.png, {pfx}_territory.png, {pfx}_misorientation.png")


def plot(M, seg, rows, bnd, fn, a, lim=None, t=None, which=("theta", "territory", "misorientation")):
    """Separate figures <base>_theta.png, <base>_territory.png, <base>_misorientation.png (fn = base + '.png', or an explicit dict
    {name: path} in `fn`).  lim = (theta_min, theta_max, misorientation_max) fixes the colour scales (animation frames: territories are
    then coloured by their mean theta, so colours mean the same in every frame); t = time in ns for the label."""
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt, matplotlib.tri as mt
    from matplotlib.collections import LineCollection, PolyCollection
    path = fn if isinstance(fn, dict) else {k: (fn[:-4] if fn.endswith(".png") else fn) + "_" + k + ".png" for k in which}
    xl, yl = (M.P[:, 0].min(), M.P[:, 0].max()), (M.P[:, 1].min(), M.P[:, 1].max())
    def fig1(title):
        f, x = plt.subplots(figsize=(6.6, 6.0), dpi=110)
        x.set_aspect("equal"); x.set_xlim(*xl); x.set_ylim(*yl); x.set_xlabel("x (um)"); x.set_ylabel("y (um)"); x.set_title(title, fontsize=10)
        if t is not None:
            x.text(0.02, 0.97, "t = %.2f ns" % t, transform=x.transAxes, va="top", fontsize=9, color="w", bbox=dict(fc="k", alpha=.5, lw=0))
        return f, x
    th0, th1, mmax = lim if lim is not None else (None, None, max([b["misorient_deg"] for b in bnd] + [1]))
    tr = mt.Triangulation(M.P[:, 0], M.P[:, 1], M.T)
    if "theta" in path:
        f, x = fig1("theta_lattice_deg"); c = x.tripcolor(tr, M.th, shading="gouraud", cmap="twilight", vmin=th0, vmax=th1)
        plt.colorbar(c, ax=x, shrink=0.8, label="deg"); f.tight_layout(); f.savefig(path["theta"]); plt.close(f)
    if "territory" in path:
        if lim is None:
            cm = plt.get_cmap("tab20"); col = np.array([cm(g % 20) for g in seg["terr"]])
        else:
            cm = plt.get_cmap("twilight"); mug = seg["mu"][np.maximum(seg["terr"], 0)]
            col = cm(np.clip((mug - th0) / (th1 - th0), 0, 1))
        col[seg["plat"] < 0] = col[seg["plat"] < 0] * 0.35 + 0.65
        f, x = fig1(f"{seg['G']} grains: plateaus (solid), ramps / cores (pale)" + ("; colour = grain mean theta" if lim is not None else ""))
        x.add_collection(PolyCollection(M.P[M.T], facecolors=col, edgecolors=col, linewidths=0.4))
        for r in rows:
            x.text(r["cx"], r["cy"], str(r["id"]), ha="center", va="center", fontsize=9, weight="bold",
                   bbox=dict(fc="w", ec="none", alpha=0.55, pad=0.6) if lim is not None else None)
        f.tight_layout(); f.savefig(path["territory"]); plt.close(f)
    if "misorientation" in path:
        f, x = fig1("boundaries between grain territories (label = misorientation, deg)")
        x.add_collection(PolyCollection(M.P[M.T], facecolors=plt.get_cmap("Greys")(0.08 + 0.12 * (seg["plat"] >= 0)), edgecolors="none"))
        segs, cc = [], []
        for b in bnd:
            for e in b["edges"]:
                segs.append(M.P[M.en[e]]); cc.append(b["misorient_deg"])
        lc = LineCollection(segs, array=np.array(cc), cmap="turbo", linewidths=2.2); lc.set_clim(0, mmax)
        x.add_collection(lc); plt.colorbar(lc, ax=x, shrink=0.8, label="misorientation (deg, raw)")
        for b in bnd:
            m = M.emid[b["edges"]].mean(0)
            x.text(m[0], m[1], f"{b['misorient_deg']:.1f}", fontsize=7, ha="center", va="center", bbox=dict(fc="w", ec="none", alpha=0.7, pad=0.5))
        f.tight_layout(); f.savefig(path["misorientation"]); plt.close(f)


if __name__ == "__main__":
    main()
