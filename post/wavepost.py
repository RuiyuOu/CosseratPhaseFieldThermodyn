"""Post-processing of dense snapshots of the finite-deformation variant-A run.
Derived fields (all documented in SMOOTHING.md):
  div v, curl v           element gradient of the P1 nodal velocity (constant per triangle)
  pressure p, shear stress  from the P0 Cauchy stress at element-mean kinematics
  GND from curl Fp        nodal recovery of the P0 Fp, optional screened-Poisson smoothing, element curl
"""
import numpy as np, meshio, scipy.sparse as sp, scipy.sparse.linalg as spl, glob, re, os, json
B_UM = 2.556e-4   # Burgers vector, um

class Mesh:
    def __init__(self, vtu):
        m = meshio.read(vtu)
        self.p = m.points[:, :2].copy(); self.tri = m.cells_dict['triangle']
        self.nn = len(self.p); self.nt = len(self.tri)
        x = self.p[self.tri, 0]; y = self.p[self.tri, 1]
        # P1 shape-function gradients: grad N_i = (b_i, c_i)/(2A)
        b = np.stack([y[:,1]-y[:,2], y[:,2]-y[:,0], y[:,0]-y[:,1]], 1)
        c = np.stack([x[:,2]-x[:,1], x[:,0]-x[:,2], x[:,1]-x[:,0]], 1)
        A2 = x[:,0]*b[:,0] + x[:,1]*b[:,1] + x[:,2]*b[:,2]
        self.A = 0.5*np.abs(A2); s = np.sign(A2)
        rows = np.repeat(np.arange(self.nt), 3); cols = self.tri.ravel()
        den = (A2).repeat(3)
        self.Gx = sp.csr_matrix(((b.ravel()/den), (rows, cols)), shape=(self.nt, self.nn))
        self.Gy = sp.csr_matrix(((c.ravel()/den), (rows, cols)), shape=(self.nt, self.nn))
        # nodal recovery R: P0 -> P1, area-weighted patch average (lumped L2 projection)
        W = sp.csr_matrix((np.repeat(self.A, 3), (cols, rows)), shape=(self.nn, self.nt))
        self.R = sp.diags(1.0/np.asarray(W.sum(1)).ravel()) @ W
        # P1 stiffness / mass for the screened-Poisson smoother
        K = sp.csr_matrix((self.nn, self.nn)); Mm = sp.csr_matrix((self.nn, self.nn))
        Kloc = np.zeros((self.nt, 3, 3)); gx = b/ A2[:,None]; gy = c/A2[:,None]
        Kloc = (gx[:,:,None]*gx[:,None,:] + gy[:,:,None]*gy[:,None,:]) * self.A[:,None,None]
        Mloc = (np.ones((3,3)) + np.eye(3))[None] * (self.A/12.0)[:,None,None]
        I = np.repeat(self.tri, 3, axis=1).ravel(); J = np.tile(self.tri, (1,3)).ravel()
        self.K = sp.csr_matrix((Kloc.ravel(), (I, J)), shape=(self.nn, self.nn))
        self.M = sp.csr_matrix((Mloc.ravel(), (I, J)), shape=(self.nn, self.nn))
        self._lu = {}
    def smooth(self, q, ell):
        """screened Poisson (I - ell^2 Laplacian) q~ = q, natural BC, P1 FEM"""
        if ell <= 0: return q
        if ell not in self._lu:
            self._lu[ell] = spl.splu((self.M + ell**2*self.K).tocsc())
        return self._lu[ell].solve(self.M @ q)

def load(vtu, mesh):
    m = meshio.read(vtu)
    ti = [i for i, c in enumerate(m.cells) if c.type == 'triangle'][0]
    g = lambda k: m.cell_data[k][ti].ravel()
    d = dict(u=m.point_data['u'][:, :2], v=m.point_data['v'][:, :2],
             thdot=m.point_data['theta_dot'].ravel(), dth=m.point_data['dtheta_deg'].ravel(),
             T=m.point_data['T'].ravel(), phi=m.point_data['phi'].ravel())
    for k in ['sig_11','sig_22','sig_12','sig_21','Fp_11','Fp_12','Fp_21','Fp_22','rhoG','slip_path_eff','slip_path_phys','rhoS','wstar','D_rate','J','p1','p2','p3','Label']:
        d[k] = g(k)
    return d

def derived(mesh, d, ell=0.04):
    o = {}
    vx, vy = d['v'][:, 0], d['v'][:, 1]
    o['divv'] = mesh.Gx @ vx + mesh.Gy @ vy           # 1/ns
    o['curlv'] = mesh.Gx @ vy - mesh.Gy @ vx          # 1/ns
    o['press'] = -0.5*(d['sig_11'] + d['sig_22'])*1e3   # MPa (sigma in GPa)
    o['tauS'] = 0.5*(d['sig_12'] + d['sig_21'])*1e3     # MPa
    o['umag'] = np.hypot(d['u'][:, 0], d['u'][:, 1])*1e3   # nm
    o['dT'] = d['T'] - 300.0
    # GND from curl Fp: (curl Fp)_{i} = d_1 Fp_{i2} - d_2 Fp_{i1}
    def rec(name): return mesh.smooth(mesh.R @ d[name], ell)
    F11, F12, F21, F22 = [rec(n) for n in ('Fp_11','Fp_12','Fp_21','Fp_22')]
    c1 = mesh.Gx @ F12 - mesh.Gy @ F11
    c2 = mesh.Gx @ F22 - mesh.Gy @ F21
    o['gnd_curl'] = np.hypot(c1, c2)/B_UM*1e12          # m^-2  (alpha magnitude / b)
    o['gnd_curl_x'], o['gnd_curl_y'] = c1, c2
    F11r, F12r, F21r, F22r = [mesh.R @ d[n] for n in ('Fp_11','Fp_12','Fp_21','Fp_22')]
    o['gnd_curl_raw'] = np.hypot(mesh.Gx @ F12r - mesh.Gy @ F11r, mesh.Gx @ F22r - mesh.Gy @ F21r)/B_UM*1e12
    o['gnd_solver'] = d['rhoG']*1e12
    for k in ('slip_path_eff','slip_path_phys','p1','p2','p3'): o[k] = d[k]
    return o

def frame_list(dirpath, base='A_dyn'):
    fs = []
    rel = os.path.join(dirpath, base+'_relaxed.vtu')
    if os.path.exists(rel): fs.append((0.0, rel))
    pvd = os.path.join(dirpath, base+'.pvd')
    if os.path.exists(pvd):
        for t, f in re.findall(r'timestep="([^"]+)" file="([^"]+)"', open(pvd).read()):
            path = os.path.join(dirpath, os.path.basename(f))
            if os.path.exists(path) and 'relaxed' not in path: fs.append((float(t), path))
    fs.sort(key=lambda a: a[0])
    # drop duplicates in time
    out = []; last = -1
    for t, f in fs:
        if abs(t - last) > 1e-6: out.append((t, f)); last = t
    return out
