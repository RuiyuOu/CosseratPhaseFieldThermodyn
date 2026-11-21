#!/usr/bin/env python3
"""B1 -- where and when does a lattice-curvature indicator ask for refinement?  2D post-processing of a run's VTU frames.

For every frame of out_<tag>/<tag>.pvd (point data theta_lattice_deg = theta_0 + theta):
  kappa_e   = |grad theta_lattice| per triangle (P1 gradient, deg/um), wrapped so that a jump across the +-90 deg branch cut is not counted
  h_target  = clamp(Delta* / kappa_e, hmin, hmax)                       (Delta* = largest orientation step allowed across one element)
  then on a raster: buffer b (every point within b of an element that wants size h gets h) and gradation g (h may grow by at most g per um:
  h(x) <= h(y) + g |x - y|), i.e. h_final(x) = min over levels h_k of h_k + g * max(0, dist(x, {h_target <= h_k}) - b)
reported per frame and per (Delta*, b):
  f_min     area fraction at h_min (h_final < 1.5 h_min)
  n2d       ideal triangle count of the graded 2D mesh  (sum dA / (sqrt(3)/4 h^2))
  n3d_10um  projected tet count for a 10 um cube with the same statistics (volume fraction = area fraction, Delesse):
            V * mean(1 / v_tet(h)), v_tet(h) = h^3 / (6 sqrt 2); uniform h_min for comparison.  Real meshers make ~20% more.
  front_um  how far from the loaded edge (x = Lx) the lattice has rotated by more than ROT deg since the first frame (99th percentile)
usage:  python3 b1_curvature.py out_<tag> [--hmin 0.1] [--hmax 0.4] [--grade 0.3] [--dstar 0.25 0.5 1] [--buffer 0 0.25 0.5 1]
        [--every 1] [--rot 1] [--pix 0.02] [--out b1_<tag>]
writes <out>/b1_frames.csv (one row per frame and Delta*, b), <out>/b1_summary.png, and <out>/b1_hfinal_<frame>.png for a few frames."""
import sys, os, re, glob, argparse, math, csv
import numpy as np, meshio
from scipy import ndimage
ap = argparse.ArgumentParser()
ap.add_argument("run"); ap.add_argument("--hmin", type=float, default=0.1); ap.add_argument("--hmax", type=float, default=0.4)
ap.add_argument("--grade", type=float, default=0.3); ap.add_argument("--dstar", type=float, nargs="+", default=[0.25, 0.5, 1.0])
ap.add_argument("--buffer", type=float, nargs="+", default=[0.0, 0.25, 0.5, 1.0]); ap.add_argument("--every", type=int, default=1)
ap.add_argument("--rot", type=float, default=1.0); ap.add_argument("--pix", type=float, default=0.02); ap.add_argument("--out", default=None)
ap.add_argument("--nlev", type=int, default=12)
a = ap.parse_args()
run = a.run.rstrip("/")
pvd = glob.glob(os.path.join(run, "*.pvd"))
frames = []
if pvd:
    for m in re.finditer(r'timestep="([^"]+)"[^>]*file="([^"]+)"', open(pvd[0]).read()):
        f = os.path.join(run, m.group(2))
        if os.path.exists(f): frames.append((float(m.group(1)), f))
if not frames: sys.exit("no frames found via %s/*.pvd" % run)
frames.sort(); frames = [fr for i, fr in enumerate(frames) if i == 0 or fr[0] != frames[i - 1][0]][::a.every]   # one frame per time
out = a.out or ("b1_" + os.path.basename(run).replace("out_", ""))
os.makedirs(out, exist_ok=True)
M0 = meshio.read(frames[0][1]); P = M0.points[:, :2]; T = M0.cells_dict["triangle"]
x0, y0 = P.min(0); x1, y1 = P.max(0); Lx, Ly = x1 - x0, y1 - y0; A_tot = Lx*Ly
# per-triangle P1 gradient operator
p0, p1, p2 = P[T[:, 0]], P[T[:, 1]], P[T[:, 2]]
det = (p1[:, 0] - p0[:, 0])*(p2[:, 1] - p0[:, 1]) - (p2[:, 0] - p0[:, 0])*(p1[:, 1] - p0[:, 1]); area = 0.5*np.abs(det)
b = np.stack([p1[:, 1] - p2[:, 1], p2[:, 1] - p0[:, 1], p0[:, 1] - p1[:, 1]], 1)/det[:, None]
c = np.stack([p2[:, 0] - p1[:, 0], p0[:, 0] - p2[:, 0], p1[:, 0] - p0[:, 0]], 1)/det[:, None]
cen = (p0 + p1 + p2)/3
hedge = np.max(np.stack([np.linalg.norm(p1 - p0, axis=1), np.linalg.norm(p2 - p1, axis=1), np.linalg.norm(p0 - p2, axis=1)], 1), 1)
def grad(v):   # v nodal; values of a triangle are unwrapped against its first node (theta lives mod 180 deg)
    w = v[T]; w = w - 180.0*np.round((w - w[:, :1])/180.0)
    return (w*b).sum(1), (w*c).sum(1)
# raster for the buffer / gradation step
nx, ny = int(math.ceil(Lx/a.pix)), int(math.ceil(Ly/a.pix))
ix = np.clip(((cen[:, 0] - x0)/a.pix).astype(int), 0, nx - 1); iy = np.clip(((cen[:, 1] - y0)/a.pix).astype(int), 0, ny - 1)
levels = np.geomspace(a.hmin, a.hmax, a.nlev)
def hfinal(ht, buf):
    grid = np.full((ny, nx), a.hmax)
    np.minimum.at(grid, (iy, ix), ht)            # smallest target in each pixel (pixels without a centroid keep hmax)
    hf = np.full((ny, nx), a.hmax)
    for hk in levels:
        mask = grid <= hk*(1 + 1e-9)
        if not mask.any(): continue
        d = ndimage.distance_transform_edt(~mask)*a.pix
        hf = np.minimum(hf, hk + a.grade*np.maximum(0.0, d - buf))
    return np.clip(hf, a.hmin, a.hmax)
vtet = lambda h: h**3/(6*math.sqrt(2))
th0 = M0.point_data["theta_lattice_deg"].ravel()
rows = []
show = set(np.linspace(0, len(frames) - 1, 4).astype(int))
for k, (t, f) in enumerate(frames):
    M = meshio.read(f); th = M.point_data["theta_lattice_deg"].ravel()
    gx, gy = grad(th); kap = np.hypot(gx, gy)                    # deg/um
    rot = th - th0; rot = rot - 180.0*np.round(rot/180.0)
    moved = np.abs(rot) > a.rot
    front = float(np.percentile(x1 - P[moved, 0], 99)) if moved.sum() > 20 else 0.0
    for ds in a.dstar:
        ht = np.clip(ds/np.maximum(kap, 1e-12), a.hmin, a.hmax)
        for buf in a.buffer:
            hf = hfinal(ht, buf)
            fmin = float((hf < 1.5*a.hmin).mean())
            n2d = float((a.pix**2/(math.sqrt(3)/4*hf**2)).sum())
            n3d = float(1000.0*np.mean(1.0/vtet(hf)))
            rows.append(dict(t_ns=t, frame=os.path.basename(f), dstar_deg=ds, buffer_um=buf, f_min=fmin, n2d=n2d, n3d_10um=n3d, front_um=front,
                             kappa_p99=float(np.percentile(kap, 99)), dmax_elem_deg=float((kap*hedge).max())))
            if k in show and ds == a.dstar[len(a.dstar)//2] and buf == a.buffer[min(2, len(a.buffer) - 1)]:
                try:
                    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
                    fig, ax = plt.subplots(1, 2, figsize=(10, 4.4))
                    im0 = ax[0].tripcolor(P[:, 0], P[:, 1], T, facecolors=np.log10(np.maximum(kap, 1e-2)), cmap="magma", vmin=-1, vmax=3)
                    ax[0].set_title("log10 |grad theta_lattice| (deg/um), t = %g ns" % t); fig.colorbar(im0, ax=ax[0])
                    im1 = ax[1].imshow(hf, origin="lower", extent=[x0, x1, y0, y1], cmap="viridis_r", vmin=a.hmin, vmax=a.hmax)
                    ax[1].set_title("target size (um): D*=%g deg, buffer %g um, f_min %.2f" % (ds, buf, fmin)); fig.colorbar(im1, ax=ax[1])
                    for x_ in ax: x_.set_aspect("equal")
                    fig.tight_layout(); fig.savefig(os.path.join(out, "b1_hfinal_%04d.png" % k), dpi=110); plt.close(fig)
                except Exception as e: print("plot failed:", e)
    print("t %6.2f ns  front %5.2f um  kappa_p99 %8.1f deg/um  " % (t, front, np.percentile(kap, 99)) +
          "  ".join("D*%g b%g: f %.2f n3d %.2fM" % (r["dstar_deg"], r["buffer_um"], r["f_min"], r["n3d_10um"]/1e6)
                    for r in rows[-len(a.dstar)*len(a.buffer):] if r["buffer_um"] == a.buffer[min(2, len(a.buffer) - 1)]), flush=True)
with open(os.path.join(out, "b1_frames.csv"), "w") as fo:
    w = csv.DictWriter(fo, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
try:
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
    for ds in a.dstar:
        for buf in a.buffer:
            rr = [r for r in rows if r["dstar_deg"] == ds and r["buffer_um"] == buf]
            tt = [r["t_ns"] for r in rr]
            ax[0].plot(tt, [r["f_min"] for r in rr], label="D*=%g, b=%g" % (ds, buf))
            ax[1].plot(tt, [r["n3d_10um"]/1e6 for r in rr])
    tt = sorted(set(r["t_ns"] for r in rows)); fr = {r["t_ns"]: r["front_um"] for r in rows}
    ax[2].plot(tt, [fr[x] for x in tt])
    ax[1].axhline(1000/vtet(a.hmin)/1e6, ls="--", c="k", lw=0.8); ax[1].text(0, 1000/vtet(a.hmin)/1e6, " uniform h_min", va="bottom")
    ax[0].set(xlabel="t (ns)", ylabel="area fraction at h_min", title="h_min %g, h_max %g, grading %g" % (a.hmin, a.hmax, a.grade))
    ax[1].set(xlabel="t (ns)", ylabel="projected tets, 10 um cube (millions)", yscale="log")
    ax[2].set(xlabel="t (ns)", ylabel="rotation front from loaded edge (um)", title="|d theta_lattice| > %g deg" % a.rot)
    ax[0].legend(fontsize=7, ncol=2); fig.tight_layout(); fig.savefig(os.path.join(out, "b1_summary.png"), dpi=110)
except Exception as e: print("summary plot failed:", e)
print("wrote", out)
