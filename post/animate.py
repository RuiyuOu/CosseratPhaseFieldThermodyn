import sys, os, numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt, matplotlib.tri as mtri
from matplotlib.colors import LogNorm
from multiprocessing import Pool
import wavepost as w
from PIL import Image
D = sys.argv[1]; OUT = sys.argv[2]; ELL = float(sys.argv[3]) if len(sys.argv) > 3 else 0.04
MAG = 400.0   # replaced by auto value in limits(): max displaced shift = 8% of block size
FR = w.frame_list(D, os.environ.get('BASE', 'A_dyn')); M = w.Mesh(FR[0][1])
# name: (title, key, cmap, symmetric?, log?, unit)
ANIMS = {
 'phi':      ('Phase field $\\phi$ (order parameter)', 'phi', 'viridis', False, False, ''),
 'pressure': ('Pressure wave: $p=-\\frac{1}{2}(\\sigma_{11}+\\sigma_{22})$', 'press', 'RdBu_r', True, False, 'MPa'),
 'shear':    ('Shear wave: $\\frac{1}{2}(\\sigma_{12}+\\sigma_{21})$', 'tauS', 'PuOr', True, False, 'MPa'),
 'divv':     ('Dilatational wave: $\\nabla\\!\\cdot v$', 'divv', 'RdBu_r', True, False, '1/ns'),
 'curlv':    ('Rotational (shear) wave: $\\nabla\\times v$', 'curlv', 'PuOr', True, False, '1/ns'),
 'microrot': ('Microrotation wave: $\\dot\\theta$', 'thdot', 'seismic', True, False, 'rad/ns'),
 'temperature': ('Temperature rise $T-300$', 'dT', 'inferno', False, False, 'K'),
 'deformation': ('Deformation $|u|$ (displacement magnified, see label)', 'umag', 'cividis', False, False, 'nm'),
 'slip_acc': ('Accumulated slip $\\sum_\\alpha\\int|\\dot\\gamma_\\alpha|dt$ (all systems)', 'slip_path_eff', 'turbo', False, True, 'shear strain'),
 'slip_sys1': ('Accumulated signed slip $\\gamma_1$ (system 1)', 'p1', 'RdBu_r', True, False, 'shear strain'),
 'slip_sys2': ('Accumulated signed slip $\\gamma_2$ (system 2)', 'p2', 'RdBu_r', True, False, 'shear strain'),
 'slip_sys3': ('Accumulated signed slip $\\gamma_3$ (system 3)', 'p3', 'RdBu_r', True, False, 'shear strain'),
 'gnd_curl': ('GND from curl $F^p$: $|\\mathrm{curl}F^p|/b$', 'gnd_curl', 'magma', False, True, 'm$^{-2}$'),
 'gnd_solver': ('GND from solver: $\\phi^2|\\nabla\\theta-\\nabla\\theta_R|/b$', 'gnd_solver', 'magma', False, True, 'm$^{-2}$'),
}
NODAL = {'phi','thdot','dT','umag'}
def get(i):
    t, f = FR[i]; d = w.load(f, M); o = w.derived(M, d, ELL); o['phi'] = d['phi']; o['thdot'] = d['thdot']
    o['u'] = d['u']; return t, o
def limits():
    idx = np.unique(np.linspace(0, len(FR)-1, min(len(FR), 12)).astype(int)); acc = {k: [] for k in ANIMS}; acc_u = []
    for i in idx:
        _, o = get(i); acc_u.append(o['u'].ravel())
        for k, (_, key, *_r) in ANIMS.items(): acc[k].append(o[key])
    global MAG
    umax = max(np.abs(np.concatenate(acc_u)).max(), 1e-12); MAG = min(0.08*10.0/umax, 2000.0)
    lim = {}
    for k, (_, key, cm, sym, lg, _u) in ANIMS.items():
        a = np.concatenate(acc[k])
        if lg: lim[k] = ((max(np.percentile(a[a > 0], 2), 1e8) if k.startswith('gnd') else max(np.percentile(a[a > 0], 5), 1e-6)) if (a > 0).any() else 1e-6, max(np.percentile(a, 99.8), 1e-5))
        elif sym: m = np.percentile(np.abs(a), 99.5); lim[k] = (-m, m)
        else: lim[k] = (np.percentile(a, 0.2), np.percentile(a, 99.8))
    return lim
LIM = None
def render(i):
    t, o = get(i); res = []
    for k, (ttl, key, cm, sym, lg, un) in ANIMS.items():
        lo, hi = LIM[k]; v = o[key]
        fig, ax = plt.subplots(figsize=(5.6, 5.0), dpi=110)
        P = M.p.copy()
        if k == 'deformation': P = P + MAG*o['u']   # u in um already? (u stored um; MAG*u)
        tr = mtri.Triangulation(P[:, 0], P[:, 1], M.tri)
        kw = dict(cmap=cm, norm=LogNorm(lo, hi)) if lg else dict(cmap=cm, vmin=lo, vmax=hi)
        pc = ax.tripcolor(tr, v, **kw, shading='flat') if len(v) == M.nt else ax.tripcolor(tr, v, **kw, shading='gouraud')
        ax.set_aspect('equal'); ax.set_xlabel('x (µm)'); ax.set_ylabel('y (µm)')
        ax.set_title(ttl, fontsize=9); ax.text(0.02, 0.97, 't = %.2f ns' % t, transform=ax.transAxes, va='top', fontsize=9, color='w', bbox=dict(fc='k', alpha=.5, lw=0))
        if k == 'deformation':
            ax.plot([0,10,10,0,0],[0,0,10,10,0],'w--',lw=.8); ax.set_xlim(-1.2,11.2); ax.set_ylim(-1.2,11.2)
            ax.set_title('Deformation $|u|$, displacement magnified %.0fx (dashed = undeformed)' % MAG, fontsize=8)
        cb = fig.colorbar(pc, ax=ax, shrink=.85); cb.set_label(un)
        fig.tight_layout(); d = os.path.join(OUT, k, 'png'); os.makedirs(d, exist_ok=True)
        fn = os.path.join(d, '%s_%04d.png' % (k, i)); fig.savefig(fn); plt.close(fig); res.append(fn)
    return res
if __name__ == '__main__':
    LIM = limits(); LIM['_mag'] = (MAG, MAG)
    import json; json.dump({k: list(map(float, v)) for k, v in LIM.items()}, open(os.path.join(OUT, 'limits.json'), 'w')) if os.path.isdir(OUT) or os.makedirs(OUT) is None else 0
    def init(l):
        global LIM, MAG; LIM = l; MAG = l['_mag'][0]
    with Pool(int(os.environ.get('NP', 4)), initializer=init, initargs=(LIM,)) as p: p.map(render, range(len(FR)))
    for k in ANIMS:
        fs = sorted(os.listdir(os.path.join(OUT, k, 'png'))); ims = [Image.open(os.path.join(OUT, k, 'png', f)).convert('P', palette=Image.ADAPTIVE) for f in fs]
        ims[0].save(os.path.join(OUT, k + '.gif'), save_all=True, append_images=ims[1:], duration=80, loop=0)
