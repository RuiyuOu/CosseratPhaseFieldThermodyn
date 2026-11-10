"""Quantify nucleation-like events: (1) new disordered regions (phi dropped where it was ~1 in the reference frame),
(2) lattice-orientation pockets (misrotation relative to the grain median, developed since the reference frame).
usage: python3 nucleus_probe.py <vtu_dir> <base> <out_prefix> [dphi_thr=0.7] [dtheta_thr_deg=2] [step=1]
Outputs <out_prefix>_nucleus.csv (one row per frame, largest regions) and <out_prefix>_nucleus_last.png.
Areas in um^2, sizes as equivalent-circle diameter in nm. Reference = first frame."""
import sys, numpy as np, scipy.sparse as sp
from scipy.sparse.csgraph import connected_components
import meshio, wavepost as w
D, BASE, OUT = sys.argv[1:4]
PHI_THR = float(sys.argv[4]) if len(sys.argv) > 4 else 0.7
TH_THR = float(sys.argv[5]) if len(sys.argv) > 5 else 2.0
STEP = int(sys.argv[6]) if len(sys.argv) > 6 else 1
FR = w.frame_list(D, BASE)[::STEP]; M = w.Mesh(FR[0][1]); tri = M.tri
def elem(a): return a[tri].mean(1)
def adj():  # element adjacency through shared edges
    e = np.sort(np.vstack([tri[:, [0, 1]], tri[:, [1, 2]], tri[:, [2, 0]]]), 1); idx = np.tile(np.arange(M.nt), 3)
    o = np.lexsort((e[:, 1], e[:, 0])); e, idx = e[o], idx[o]; same = (e[1:] == e[:-1]).all(1)
    a, b = idx[:-1][same], idx[1:][same]; return sp.coo_matrix((np.ones(len(a)), (a, b)), shape=(M.nt, M.nt)).tocsr()
A = adj()
def regions(mask):
    if not mask.any(): return []
    sub = np.where(mask)[0]; n, lab = connected_components(A[sub][:, sub], directed=False); out = []
    for c in range(n):
        ids = sub[lab == c]; ar = M.A[ids].sum(); out.append((ar, ids))
    return sorted(out, key=lambda t: -t[0])
d0 = w.load(FR[0][1], M); phi0 = elem(d0['phi']); th0 = elem(d0['dth']) if False else None
th_ref = elem(np.asarray(meshio.read(FR[0][1]).point_data['theta_lattice_deg']).ravel()); grain = d0['Label'].astype(int)
rows = []; last = None
for t, f in FR:
    m = meshio.read(f); d = w.load(f, M); phi = elem(d['phi']); th = elem(np.asarray(m.point_data['theta_lattice_deg']).ravel())
    dth = th - th_ref; dev = np.zeros(M.nt)
    for g in np.unique(grain):
        s = (grain == g) & (phi > 0.95)
        if s.sum() > 20: dev[grain == g] = dth[grain == g] - np.median(dth[s])
    drop = (phi < PHI_THR) & (phi0 >= 0.95)
    pock = (np.abs(dev) > TH_THR) & (phi > 0.95)
    Tn = elem(d['T']); rho = d['rhoS']
    row = dict(t_ns=t, phi_drop_area=M.A[drop].sum(), n_phi_drop=len(regions(drop)), pocket_area=M.A[pock].sum(), max_abs_dev_deg=np.abs(dev).max(), max_T=Tn.max())
    rg = regions(drop); rp = regions(pock)
    for nm, r in (('drop', rg), ('pocket', rp)):
        if r:
            ar, ids = r[0]; c = (M.p[tri[ids]].mean(1) * M.A[ids][:, None]).sum(0) / M.A[ids].sum()
            row.update({nm+'_largest_area': ar, nm+'_largest_diam_nm': 2e3*np.sqrt(ar/np.pi), nm+'_x': c[0], nm+'_y': c[1], nm+'_meanT': Tn[ids].mean(), nm+'_mean_rhoS': rho[ids].mean()})
        else: row.update({nm+'_largest_area': 0, nm+'_largest_diam_nm': 0, nm+'_x': np.nan, nm+'_y': np.nan, nm+'_meanT': np.nan, nm+'_mean_rhoS': np.nan})
    rows.append(row); last = (t, phi, dev, drop, pock)
import csv
keys = list(rows[0].keys()); 
with open(OUT + '_nucleus.csv', 'w', newline='') as fh:
    wr = csv.DictWriter(fh, keys); wr.writeheader(); [wr.writerow(r) for r in rows]
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt, matplotlib.tri as mtri
t, phi, dev, drop, pock = last; tr = mtri.Triangulation(M.p[:, 0], M.p[:, 1], tri)
fig, ax = plt.subplots(1, 3, figsize=(15, 4.8)); fig.suptitle('%s  t=%.2f ns' % (BASE, t))
for a, v, ttl, cm, kw in ((ax[0], phi, 'phi', 'viridis', dict(vmin=0, vmax=1)), (ax[1], dev, 'lattice misrotation vs grain median since ref (deg)', 'RdBu_r', dict(vmin=-10, vmax=10)), (ax[2], 1.0*drop + 2.0*pock, 'new phi-drop (1) / orientation pocket (2)', 'cool', dict(vmin=0, vmax=3))):
    pc = a.tripcolor(tr, v, cmap=cm, shading='flat', **kw); a.set_aspect('equal'); a.set_title(ttl, fontsize=9); fig.colorbar(pc, ax=a, shrink=.8)
fig.tight_layout(); fig.savefig(OUT + '_nucleus_last.png', dpi=110)
print('wrote', OUT + '_nucleus.csv', len(rows), 'frames; last phi-drop area %.4g um2, pocket area %.4g um2' % (rows[-1]['phi_drop_area'], rows[-1]['pocket_area']))
