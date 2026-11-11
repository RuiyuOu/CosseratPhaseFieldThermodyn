"""Post-processing helpers for the solver snapshots (mesh, frame list, derived fields; documented in DATA_FIELDS.md):
  div v, curl v           element gradient of the P1 nodal velocity (constant per triangle)
  pressure p, shear stress  from the P0 Cauchy stress at element-mean kinematics
  GND density             the solver's Cermelli-Gurtin rhoG = |G - G_R|/b, G = Fp Curl Fp (cell field, converted to m^-2)
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

def load(vtu, mesh):
    m = meshio.read(vtu)
    ti = [i for i, c in enumerate(m.cells) if c.type == 'triangle'][0]
    g = lambda k: m.cell_data[k][ti].ravel()
    d = dict(u=m.point_data['u'][:, :2], v=m.point_data['v'][:, :2],
             thdot=m.point_data['theta_dot'].ravel(), dth=m.point_data['dtheta_deg'].ravel(),
             T=m.point_data['T'].ravel(), phi=m.point_data['phi'].ravel(), thl=m.point_data['theta_lattice_deg'].ravel())
    for k in ['sig_33','sig_11','sig_22','sig_12','sig_21','Fp_11','Fp_12','Fp_21','Fp_22','rhoG','slip_path_eff','slip_path_phys','rhoS','wstar','D_rate','J','p1','p2','p3','Label']:
        d[k] = g(k)
    for k in ['rho_mob', 'eta_drag']:   # drag density rho_m and low-stress drag eta_d (cell data)
        if k in m.cell_data: d[k] = g(k)
    return d

def derived(mesh, d):
    o = {}
    vx, vy = d['v'][:, 0], d['v'][:, 1]
    o['divv'] = mesh.Gx @ vx + mesh.Gy @ vy           # 1/ns
    o['curlv'] = mesh.Gx @ vy - mesh.Gy @ vx          # 1/ns
    o['press'] = -0.5*(d['sig_11'] + d['sig_22'])*1e3   # MPa (sigma in GPa)
    o['tauS'] = 0.5*(d['sig_12'] + d['sig_21'])*1e3     # MPa
    o['sigF'] = np.sqrt(d['sig_11']**2 + d['sig_22']**2 + d['sig_12']**2 + d['sig_21']**2 + d['sig_33']**2)*1e3   # Frobenius norm, MPa
    o['umag'] = np.hypot(d['u'][:, 0], d['u'][:, 1])*1e3   # nm
    o['dT'] = d['T'] - 300.0
    o['gnd'] = d['rhoG']*1e12                          # solver's Cermelli-Gurtin GND density, um^-2 -> m^-2
    for k in ('slip_path_eff','slip_path_phys','p1','p2','p3'): o[k] = d[k]
    o['thl'] = d['thl']
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
