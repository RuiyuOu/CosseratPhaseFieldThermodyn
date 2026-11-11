import sys, os, numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt, matplotlib.tri as mtri
from matplotlib.collections import LineCollection
from matplotlib.colors import LogNorm
from multiprocessing import Pool
import wavepost as w
from PIL import Image
D = sys.argv[1]; OUT = sys.argv[2]
MAG = 400.0   # replaced by auto value in limits(): max displaced shift = 8% of block size
FR = w.frame_list(D, os.environ.get('BASE', 'A_dyn')); M = w.Mesh(FR[0][1])
# name: (title, key, cmap, symmetric?, log?, unit)
ANIMS = {
 'phi':      ('Phase field $\\phi$ (order parameter)', 'phi', 'viridis', False, False, ''),
 'pressure': ('Pressure wave: $p=-\\frac{1}{2}(\\sigma_{11}+\\sigma_{22})$', 'press', 'RdBu_r', True, False, 'MPa'),
 'shear':    ('Shear wave: $\\frac{1}{2}(\\sigma_{12}+\\sigma_{21})$', 'tauS', 'PuOr', True, False, 'MPa'),
 'stress_frob': ('Local stress $\\|\\sigma\\|_F$ (Frobenius norm)', 'sigF', 'magma', False, True, 'MPa'),
 'divv':     ('Dilatational wave: $\\nabla\\!\\cdot v$', 'divv', 'RdBu_r', True, False, '1/ns'),
 'curlv':    ('Rotational (shear) wave: $\\nabla\\times v$', 'curlv', 'PuOr', True, False, '1/ns'),
 'microrot': ('Microrotation wave: $\\dot\\theta$', 'thdot', 'seismic', True, False, 'rad/ns'),
 'temperature': ('Temperature rise $T-300$', 'dT', 'inferno', False, False, 'K'),
 'deformation': ('Deformation $|u|$ (displacement magnified, see label)', 'umag', 'cividis', False, False, 'nm'),
 'slip_acc': ('Accumulated slip $\\sum_\\alpha\\int|\\dot\\gamma_\\alpha|dt$ (all systems)', 'slip_path_eff', 'turbo', False, True, 'shear strain'),
 'slip_sys1': ('Accumulated signed slip $\\gamma_1$ (system 1)', 'p1', 'RdBu_r', True, False, 'shear strain'),
 'slip_sys2': ('Accumulated signed slip $\\gamma_2$ (system 2)', 'p2', 'RdBu_r', True, False, 'shear strain'),
 'slip_sys3': ('Accumulated signed slip $\\gamma_3$ (system 3)', 'p3', 'RdBu_r', True, False, 'shear strain'),
 'rotation': ('Lattice rotation since the relaxed state: $\\theta_{lat}(t)-\\theta_{lat}(0)$', 'rot', 'RdBu_r', True, False, 'deg'),
 'lattice_slip': ('Lattice orientation + slip directions', 'thl', 'twilight_shifted', False, False, 'deg'),
 'gnd':      ('GND density (Cermelli-Gurtin): $|G-G_R|/b$, $G=F^p\\,\\mathrm{Curl}\\,F^p$', 'gnd', 'magma', False, True, 'm$^{-2}$'),
}
_d0 = w.load(FR[0][1], M); THL0 = _d0['thl'].copy()   # reference orientation for the 'rotation' field
if 'rho_mob' in _d0:   # drag fields (always written by cosserat_kwc.edp)
    ANIMS['rho_mob'] = ('Mobile dislocation density $\\rho_m$ used in the drag', 'rho_mob', 'magma', False, True, 'm$^{-2}$')
    ANIMS['eta_drag'] = ('Drag $\\eta_d=B(T)/(b^2\\rho_m)$ (low-stress, before the $v_c$ factor)', 'eta_drag', 'viridis', False, True, 'Pa s')
NODAL = {'phi','thdot','dT','umag'}
def get(i):
    t, f = FR[i]; d = w.load(f, M); o = w.derived(M, d); o['phi'] = d['phi']; o['thdot'] = d['thdot']
    o['u'] = d['u']; o['rot'] = ((o['thl'] - THL0 + 180.0) % 360.0) - 180.0
    if 'rho_mob' in d: o['rho_mob'] = d['rho_mob']*1e12; o['eta_drag'] = d['eta_drag']   # um^-2 -> m^-2; GPa ns = Pa s
    return t, o
def limits(keys=None):
    # keys: only these fields (a field added after an earlier call); None = everything incl. the displacement magnification and slip scale
    keys = list(ANIMS) if keys is None else keys; full = len(keys) == len(ANIMS)
    idx = np.unique(np.linspace(0, len(FR)-1, min(len(FR), 12)).astype(int)); acc = {k: [] for k in keys}; acc_u = []; acc_g = []
    for i in idx:
        _, o = get(i); acc_u.append(o['u'].ravel()); acc_g.append(np.abs(np.concatenate([o['p1'], o['p2'], o['p3']])))
        for k in keys: acc[k].append(o[ANIMS[k][1]])
    global MAG
    lim = {}
    if full:
        umax = max(np.abs(np.concatenate(acc_u)).max(), 1e-12); MAG = min(0.08*10.0/umax, 2000.0)
        g = np.concatenate(acc_g); lim['_g99'] = (max(np.percentile(g, 99.5), 1e-4),) * 2
    for k in keys:
        _, key, cm, sym, lg, _u = ANIMS[k]
        a = np.concatenate(acc[k])
        if k == 'lattice_slip': lim[k] = (-90.0, 90.0); continue
        if lg: lim[k] = ((max(np.percentile(a[a > 0], 2), 1e8) if k.startswith('gnd') else max(np.percentile(a[a > 0], 5), 1e-6)) if (a > 0).any() else 1e-6, max(np.percentile(a, 99.8), 1e-5))
        elif sym: m = np.percentile(np.abs(a), 99.5); m = max(m, 1.0) if k == 'rotation' else m; m = float(os.environ.get('ROTLIM', m)) if k == 'rotation' else m; lim[k] = (-m, m)
        else: lim[k] = (np.percentile(a, 0.2), np.percentile(a, 99.8))
    return lim
LIM = None
_SL = {}
ALPHA = np.degrees(np.arctan(np.sqrt(2.0)))   # lattice slip directions: 0, +54.74, -54.74 deg from the lattice x-axis (solver: slx, sly)
def slip_overlay(ax, o, ell=0.2, h=0.26):
    # regular grid -> element lookup; draw each system's direction (rotated by the local lattice angle) as a short line
    if 'tf' not in _SL:
        tr = mtri.Triangulation(M.p[:, 0], M.p[:, 1], M.tri); xs = np.arange(h/2, 10, h); X, Y = np.meshgrid(xs, xs)
        _SL.update(tf=tr.get_trifinder(), X=X.ravel(), Y=Y.ravel())
    el = _SL['tf'](_SL['X'], _SL['Y']); ok = el >= 0
    X, Y, el = _SL['X'][ok], _SL['Y'][ok], el[ok]
    th = np.radians(o['thl'][M.tri[el]].mean(1)); ph = o['phi'][M.tri[el]].mean(1); keep = ph > 0.9    # skip grain-boundary cores
    g99 = LIM['_g99'][0]; cols = ('#d62728', '#1f77b4', '#2ca02c'); segs = []
    for si, (off, key) in enumerate(((0.0, 'p1'), (ALPHA, 'p2'), (-ALPHA, 'p3'))):
        ang = th + np.radians(off); dx, dy = 0.5*ell*np.cos(ang), 0.5*ell*np.sin(ang)
        act = np.minimum(np.abs(o[key][el]) / g99, 1.0)
        P0 = np.stack([X - dx, Y - dy], 1); P1 = np.stack([X + dx, Y + dy], 1)
        lc = LineCollection(np.stack([P0, P1], 1)[keep], colors=cols[si], linewidths=(0.5 + 2.2*act)[keep], alpha=0.9, capstyle='round'); ax.add_collection(lc)
    from matplotlib.lines import Line2D
    ax.legend(handles=[Line2D([0],[0],color=c,lw=2,label=l) for c,l in zip(cols,('system 1 (0 deg)','system 2 (+54.7)','system 3 (-54.7)'))], loc='lower left', fontsize=7, framealpha=0.8)
def render(a):
    i, ks = a if isinstance(a, tuple) else (a, list(ANIMS))
    t, o = get(i); res = []
    for k in ks:
        ttl, key, cm, sym, lg, un = ANIMS[k]
        lo, hi = LIM[k]; v = o[key]
        fig, ax = plt.subplots(figsize=(5.6, 5.0), dpi=110)
        P = M.p.copy()
        if k == 'deformation': P = P + MAG*o['u']   # u in um, magnified MAG times
        tr = mtri.Triangulation(P[:, 0], P[:, 1], M.tri)
        v = ((v + 90.0) % 180.0) - 90.0 if k == 'lattice_slip' else v
        if lg: v = np.maximum(v, lo)   # zero (e.g. no GND) shows as the lowest colour, not as blank
        kw = dict(cmap=cm, norm=LogNorm(lo, hi)) if lg else dict(cmap=cm, vmin=lo, vmax=hi)
        pc = ax.tripcolor(tr, v, **kw, shading='flat') if len(v) == M.nt else ax.tripcolor(tr, v, **kw, shading='gouraud')
        ax.set_aspect('equal'); ax.set_xlabel('x (µm)'); ax.set_ylabel('y (µm)')
        ax.set_title(ttl, fontsize=9); ax.text(0.02, 0.97, 't = %.2f ns' % t, transform=ax.transAxes, va='top', fontsize=9, color='w', bbox=dict(fc='k', alpha=.5, lw=0))
        if k == 'lattice_slip': slip_overlay(ax, o)
        if k == 'deformation':
            ax.plot([0,10,10,0,0],[0,0,10,10,0],'w--',lw=.8); ax.set_xlim(-1.2,11.2); ax.set_ylim(-1.2,11.2)
            ax.set_title('Deformation $|u|$, displacement magnified %.0fx (dashed = undeformed)' % MAG, fontsize=8)
        cb = fig.colorbar(pc, ax=ax, shrink=.85); cb.set_label(un)
        fig.tight_layout(); d = os.path.join(OUT, k, 'png'); os.makedirs(d, exist_ok=True)
        fn = os.path.join(d, '%s_%04d.png' % (k, i)); fig.savefig(fn); plt.close(fig); res.append(fn)
    return res
def sig(i):
    t, f = FR[i]; return [round(t, 6), os.path.basename(f), os.path.getsize(f)]
def png(k, i): return os.path.join(OUT, k, 'png', '%s_%04d.png' % (k, i))
if __name__ == '__main__':
    # Incremental: frames whose PNGs already exist for every field (same time/file/size as in manifest.json) are kept; the colour limits of
    # the first call are kept too (limits.json), so old and new frames share the scales.  A field added later (e.g. 'rotation') gets its own
    # limits and is rendered for all frames without touching the others.  FORCE=1: redo all; RELIM=1: new limits + redo all.
    import json
    os.makedirs(OUT, exist_ok=True); LF = os.path.join(OUT, 'limits.json'); MF = os.path.join(OUT, 'manifest.json')
    force = os.environ.get('FORCE', '0') == '1' or os.environ.get('RELIM', '0') == '1'
    LIM = None
    if os.path.exists(LF) and not force:
        LIM = {k: tuple(v) for k, v in json.load(open(LF)).items()}
        if '_mag' not in LIM or '_g99' not in LIM or not (set(ANIMS) & set(LIM)): LIM = None
    fresh = LIM is None
    if fresh:
        LIM = limits(); LIM['_mag'] = (MAG, MAG); json.dump({k: list(map(float, v)) for k, v in LIM.items()}, open(LF, 'w'))
    else:
        MAG = LIM['_mag'][0]
        newk = [k for k in ANIMS if k not in LIM]
        if newk:
            print('animate: new field(s) %s: computing their colour limits' % ', '.join(newk), flush=True)
            LIM.update(limits(newk)); json.dump({k: list(map(float, v)) for k, v in LIM.items()}, open(LF, 'w'))
    old = {} if (fresh or not os.path.exists(MF)) else json.load(open(MF))
    todo = []
    for i in range(len(FR)):
        ks = list(ANIMS) if old.get(str(i)) != sig(i) else [k for k in ANIMS if not os.path.exists(png(k, i))]
        if ks: todo.append((i, ks))
    nfr = len({i for i, _ in todo})
    print('animate: %d frames, %d already complete, %d to render%s' % (len(FR), len(FR) - nfr, nfr, ' (new colour limits)' if fresh else ''), flush=True)
    def init(l):
        global LIM, MAG; LIM = l; MAG = l['_mag'][0]
    if todo:
        with Pool(int(os.environ.get('NP', 4)), initializer=init, initargs=(LIM,)) as p: p.map(render, todo)
    json.dump({str(i): sig(i) for i in range(len(FR))}, open(MF, 'w'))
    touched = set(k for _, ks in todo for k in ks)
    for k in ANIMS:
        gf = os.path.join(OUT, k + '.gif')
        if k not in touched and os.path.exists(gf): continue                # nothing new for this field: keep its GIF
        fs = ['%s_%04d.png' % (k, i) for i in range(len(FR))]
        ims = [Image.open(os.path.join(OUT, k, 'png', f)).convert('P', palette=Image.ADAPTIVE) for f in fs]
        ims[0].save(gf, save_all=True, append_images=ims[1:], duration=80, loop=0)
