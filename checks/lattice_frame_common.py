"""Shared pointwise model for checks/test_stress_follows_lattice.py and checks/test_slip_follows_lattice.py.

Reproduces the production constitutive kinematics at one material point (2D plane strain, numpy only):
  F = R Ubar Fp,  R = polar(F Fp^-1) (exact constraint),  multiplier lambda from the torque balance,
  P = R (Pi + 2 lambda da/dUbar) Fp^-T,  tau = P F^T,  M = Ubar^T Pi,  tau_a = M : (s_a x n_a).
Elastic laws: 'svk'  psi = 1/2 E:C:E + 2 mu_c a^2,  E = (Ubar^T Ubar - I)/2  (plain St Venant-Kirchhoff, default)
              'prod' psi = Vinet W(J) + 1/2 Et:C:Et + 2 mu_c a^2,  Et = (J^-2/3 Ubar^T Ubar - I)/2  (production law)
Crystal data: cubic Cu in the solver's lattice frame x=[110], y=[001], z=[-110]; slip directions as in the solver
(s = (1,0), (cos aR, sin aR), (cos aR, -sin aR), aR = atan(sqrt 2), n = (-s_y, s_x)).
Derivatives by complex step (exact to round-off).  Python >= 3.8, numpy.
"""
import os, re
import numpy as np

C11, C12, C44 = 168.4, 121.4, 75.4           # Cu, cube axes (GPa)
MUC = 754.0                                  # production mu_c (GPa)
K0, KP = (C11 + 2*C12)/3, 5.0                # Vinet (production defaults)
AR = np.arctan(np.sqrt(2.0))
SLIP_C = [np.array([1.0, 0.0]), np.array([np.cos(AR), np.sin(AR)]), np.array([np.cos(AR), -np.sin(AR)])]
I3 = np.eye(3)
WHAT = np.array([[0.0, -1.0], [1.0, 0.0]])
DADU = 0.5*WHAT                              # d a / d Ubar,  a = (Ubar21 - Ubar12)/2


def star(Q, C):
    return np.einsum('ia,jb,kc,ld,abcd->ijkl', Q, Q, Q, Q, C)


def _cubic():
    C = np.zeros((3, 3, 3, 3))
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    C[i, j, k, l] = (C12*(i == j)*(k == l) + C44*((i == k)*(j == l) + (i == l)*(j == k))
                                     + (C11 - C12 - 2*C44)*(i == j == k == l))
    return C


_LF = np.array([[1, 1, 0], [0, 0, 1], [-1, 1, 0]], float)
_LF /= np.linalg.norm(_LF, axis=1)[:, None]
CC = star(_LF, _cubic())                     # crystal tensor C^c in the solver's lattice frame


def R2(t):
    c, s = np.cos(t), np.sin(t)
    return np.array([[c, -s], [s, c]])


def R3(t):
    Q = np.eye(3); Q[:2, :2] = R2(t); return Q


def rot_angle(A):
    """exact 2D polar angle of A (det A > 0)"""
    return np.arctan2(A[1, 0] - A[0, 1], A[0, 0] + A[1, 1])


def slip_systems(theta0):
    """(s_a, n_a) of the three systems for a grain of orientation theta0 (rotated once, as the solver does at setup)"""
    Q = R2(theta0)
    return [(Q @ s, Q @ np.array([-s[1], s[0]])) for s in SLIP_C]


def slip_shear(gamma, a, theta0=0.0):
    s, n = slip_systems(theta0)[a]
    return np.eye(2) + gamma*np.outer(s, n)


def psi(U, C, law="svk"):
    U3 = np.zeros((3, 3), dtype=U.dtype); U3[:2, :2] = U; U3[2, 2] = 1
    Cb = U3.T @ U3
    a = 0.5*(U[1, 0] - U[0, 1])
    if law == "svk":
        E = 0.5*(Cb - I3)
        return 0.5*np.einsum('ij,ijkl,kl', E, C, E) + 2*MUC*a**2
    J = U[0, 0]*U[1, 1] - U[0, 1]*U[1, 0]
    Et = 0.5*(np.exp(-2.0/3.0*np.log(J))*Cb - I3)
    eta = 1.5*(KP - 1); x = np.exp(np.log(J)/3)
    Wv = 9*K0/eta**2*(1 - (1 - eta*(1 - x))*np.exp(eta*(1 - x)))
    return 0.5*np.einsum('ij,ijkl,kl', Et, C, Et) + Wv + 2*MUC*a**2


def dpsi(U, C, law="svk", h=1e-30):
    G = np.zeros((2, 2))
    for i in range(2):
        for j in range(2):
            E = np.zeros((2, 2), complex); E[i, j] = 1j*h
            G[i, j] = psi(U.astype(complex) + E, C, law).imag/h
    return G


def state(F, Fp, C, law="svk"):
    """production stress state at (F, Fp) with elasticity tensor C (reference components). Returns dict with
    tau = P F^T, M (Mandel), Fe, R, Ubar, lambda."""
    Fe = F @ np.linalg.inv(Fp)
    R = R2(rot_angle(Fe)); U = R.T @ Fe
    Pi = dpsi(U, C, law)
    dU = -WHAT @ U                                            # d Ubar / d theta
    lam = -np.sum(Pi*dU)/(2*np.sum(DADU*dU))                  # torque balance d/dtheta [psi + 2 lam a] = 0
    Pit = Pi + 2*lam*DADU
    P = R @ Pit @ np.linalg.inv(Fp).T
    return {"tau": P @ F.T, "M": U.T @ Pit, "Fe": Fe, "R": R, "U": U, "lam": lam}


def resolved(st, systems):
    """tau_a = M : (s_a x n_a) with the given (fixed) slip vectors, plus the current lattice vectors Fe s, Fe^-T n"""
    Fe = st["Fe"]; Fit = np.linalg.inv(Fe).T
    return [(s @ st["M"] @ n, Fe @ s, Fit @ n) for s, n in systems]


def load_set(seed=0, nrand=20):
    """35 stretch states: uniaxial x / y and equibiaxial at 0.02, 0.1, 0.25, -0.2; pure shears; random SPD (J ~ 0.75-1.3)"""
    rng = np.random.default_rng(seed); V = []
    for e in (0.02, 0.1, 0.25, -0.2):
        V += [np.diag([1 + e, 1]), np.diag([1, 1 + e]), np.diag([1 + e, 1 + e])]
    for g in (0.05, 0.2, 0.4):
        t = np.arctanh(g); V.append(np.array([[np.cosh(t), np.sinh(t)], [np.sinh(t), np.cosh(t)]]))
    for _ in range(nrand):
        A = rng.normal(size=(2, 2))*0.25; A = 0.5*(A + A.T); w, Q = np.linalg.eigh(A)
        V.append(Q @ np.diag(np.exp(w)) @ Q.T)
    return V


ALPHAS = np.deg2rad([0, 10, 22.5, 30, 45, 60, 77, 90, 135, -40, -100, 170])


def history_Fp(alpha):
    """T3 plastic history: slips on the three solver systems (gamma 0.5, -0.2, 0.3) followed by a pure spin R(-alpha)"""
    return slip_shear(0.3, 1) @ slip_shear(-0.2, 2) @ slip_shear(0.5, 0) @ R2(-alpha)


def stepped_spin(V, wtot_deg=-73.0, nst=200, every=20):
    """T4: F(t) log-linear from I to V; Fp <- exp(dw* What) Fp each step (exact, as trialFp). Yields (F, Fp, alpha)."""
    w, Q = np.linalg.eigh(V); Fp = np.eye(2); dw = np.deg2rad(wtot_deg)/nst
    for n in range(1, nst + 1):
        Fp = (np.cos(dw)*np.eye(2) + np.sin(dw)*WHAT) @ Fp
        if n % every == 0:
            yield Q @ np.diag(np.exp(np.log(w)*n/nst)) @ Q.T, Fp.copy(), -dw*n


def check_against_solver(edp=None):
    """compare the constants and slip directions with the solver source if it is found; returns a list of messages"""
    here = os.path.dirname(os.path.abspath(__file__))
    cands = [edp] if edp else [os.path.join(here, "..", f) for f in ("cosserat_kwc_mpi.edp", "cosserat_kwc.edp")]
    for p in cands:
        if p and os.path.isfile(p):
            src = open(p).read(); msgs = []
            want = {"Clxx": CC[0, 0, 0, 0], "Clyy": CC[1, 1, 1, 1], "Clxy": CC[0, 0, 1, 1], "Clss": CC[0, 1, 0, 1],
                    "Clxz": CC[0, 0, 2, 2], "Clyz": CC[1, 1, 2, 2]}
            for k, v in want.items():
                m = re.search(r"\b" + k + r"\s*=\s*([-\d.eE+]+)", src)
                if m is None:
                    msgs.append(("WARN", f"{k} not found in {os.path.basename(p)}"))
                elif abs(float(m.group(1)) - v) > 0.06:
                    msgs.append(("FAIL", f"{k} = {m.group(1)} in the solver, {v:.2f} here"))
                else:
                    msgs.append(("ok", f"{k} = {m.group(1)} (here {v:.2f})"))
            if re.search(r"aR\s*=\s*atan\(sqrt\(2\.?\)\)", src) and "sly = [0., sin(aR), -sin(aR)]" in src:
                msgs.append(("ok", "slip directions (1,0), (cos aR, +-sin aR), aR = atan(sqrt 2), n = (-s_y, s_x)"))
            else:
                msgs.append(("WARN", "slip-direction definition not recognised in the solver source"))
            return p, msgs
    return None, [("WARN", "solver source not found; constants not compared")]
