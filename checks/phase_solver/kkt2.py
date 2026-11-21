#!/usr/bin/env python3
# Phase-field step of cosserat_kwc.edp with explicit bound multipliers (lambda_lo for phi >= lo, lambda_hi for phi <= hi):
#   r(phi_k) + A (phi - phi_k) - lambda_lo + lambda_hi = 0,  lambda_lo >= 0, phi >= lo, lambda_lo (phi - lo) = 0,  lambda_hi >= 0, phi <= hi, lambda_hi (hi - phi) = 0
# solved by the primal-dual active-set method (Hintermueller, Ito, Kunisch 2002) with different linear solvers for the KKT system,
# against the solver's own approach (active nodes eliminated, reduced SPD system).  Same 10-step sequence and production mesh as bench.edp.
import sys, time, numpy as np, scipy.sparse as sp, scipy.sparse.linalg as sla, qdldl
# ---------------- mesh, smoothed angle, constants ----------------
with open("mesh.msh") as f:
    nv, nt, nb = map(int, f.readline().split())
    V = np.array([f.readline().split() for _ in range(nv)], float)[:, :2]
    Tr = np.array([f.readline().split() for _ in range(nt)], int)
tri, lab = Tr[:, :3] - 1, Tr[:, 3]
thg = np.loadtxt("grain_angles_rad.txt")
f0, aK, sK, dlt, dt, T0, cs = 0.005944, 0.10, 1.05, 5e-4, 0.01, 300., 0.5*48.*(2.556e-4)**2
x0, x1, x2 = V[tri[:, 0]], V[tri[:, 1]], V[tri[:, 2]]
det = (x1[:, 0]-x0[:, 0])*(x2[:, 1]-x0[:, 1]) - (x2[:, 0]-x0[:, 0])*(x1[:, 1]-x0[:, 1]); area = 0.5*np.abs(det)
gx = np.stack([(x1[:, 1]-x2[:, 1]), (x2[:, 1]-x0[:, 1]), (x0[:, 1]-x1[:, 1])], 1)/det[:, None]
gy = np.stack([(x2[:, 0]-x1[:, 0]), (x0[:, 0]-x2[:, 0]), (x1[:, 0]-x0[:, 0])], 1)/det[:, None]
s15 = np.sqrt(15.); a1, b1, a2, b2 = (6-s15)/21, (9+2*s15)/21, (6+s15)/21, (9-2*s15)/21
L = np.array([[1/3, 1/3, 1/3], [a1, a1, b1], [a1, b1, a1], [b1, a1, a1], [a2, a2, b2], [a2, b2, a2], [b2, a2, a2]])
W = np.array([9/40] + [(155-s15)/1200]*3 + [(155+s15)/1200]*3)
# fixed CSR pattern and the position of each element entry in it
rows = np.repeat(tri, 3, axis=1).ravel(); cols = np.tile(tri, (1, 3)).ravel()
P = sp.csr_matrix((np.arange(len(rows)) + 1.0, (rows, cols)), shape=(nv, nv)); P.sum_duplicates()
pat = sp.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(nv, nv)); pat.sum_duplicates(); pat.sort_indices()
key = rows.astype(np.int64)*nv + cols; pk = np.repeat(np.arange(nv), np.diff(pat.indptr)).astype(np.int64)*nv + pat.indices
pos = np.searchsorted(pk, key)
Ke = (gx[:, :, None]*gx[:, None, :] + gy[:, :, None]*gy[:, None, :])*area[:, None, None]
Me = np.einsum("q,qa,qb->ab", W, L, L)[None]*area[:, None, None]
def assemble(vals):
    v = np.bincount(pos, weights=vals.ravel(), minlength=pat.nnz); A = pat.copy(); A.data = v; return A
Mm, Km = assemble(Me), assemble(Ke)
# smoothed initial angle (eq. theta0) and N per element
Bp = sp.csr_matrix((np.repeat(area/3, 3), (tri.ravel(), np.repeat(np.arange(nt), 3))), shape=(nv, nt))
th0 = sla.spsolve((Mm + 0.04**2*Km).tocsc(), Bp @ thg[lab - 1])
gtx = (th0[tri]*gx).sum(1); gty = (th0[tri]*gy).sum(1); N = np.sqrt(gtx**2 + gty**2 + dlt**2) - dlt
kK = f0*aK**2; Kf = kK*Km
lump = np.asarray(Mm.sum(1)).ravel()
def phase_AR(phi, phik, phin, eta, rho):
    """matrix at phik and residual at phi (vphA, vphR; T = T0, so the Arrhenius factor is 1)"""
    arr = eta/dt
    pq = phi[tri] @ L.T; kq = phik[tri] @ L.T; nq = phin[tri] @ L.T       # nt x 7
    cq = arr + f0 + f0*sK*N[:, None]*2/(1 - kq)**2
    fq = arr*(pq - nq) + f0*pq + f0*sK*N[:, None]*2*pq/(1 - pq) - (f0 - cs*rho)
    wA = (cq*W)*area[:, None]; wR = (fq*W)*area[:, None]
    Me_ = np.einsum("eq,qa,qb->eab", wA, L, L)
    A = assemble(Me_) + Kf
    R = np.bincount(tri.ravel(), weights=(wR @ L).ravel(), minlength=nv) + Kf @ phi
    return A, R
# ---------------- linear solvers ----------------
def cg_jacobi(A, b, x0=None):
    it = [0]; y, info = sla.cg(A, b, x0=x0, rtol=1e-13, M=sp.diags(1/A.diagonal()), callback=lambda xk: it.__setitem__(0, it[0]+1)); return y, it[0]
def solve_elim(A, rhs, lo_set, hi_set, lo, hi, solver, phik):
    act = lo_set | hi_set; fr = np.flatnonzero(~act); ia = np.flatnonzero(act)
    x = np.where(hi_set, hi, lo).astype(float); x[fr] = 0
    Aff = A[fr][:, fr]; bf = rhs[fr] - A[fr][:, ia] @ x[ia]
    if solver == "elim-cg": y, it = cg_jacobi(Aff.tocsr(), bf, phik[fr])
    else: y, it = sla.splu(Aff.tocsc()).solve(bf), 0
    x[fr] = y; lam = A @ x - rhs                                     # multipliers recovered: lambda_lo - lambda_hi = A x - rhs (= -r_lin)
    return x, np.where(lo_set, lam, 0.), np.where(hi_set, -lam, 0.), it
def solve_kkt(A, rhs, lo_set, hi_set, lo, hi, solver, phik):
    il, ih = np.flatnonzero(lo_set), np.flatnonzero(hi_set); nl, nh = len(il), len(ih)
    El = sp.csr_matrix((np.ones(nl), (np.arange(nl), il)), shape=(nl, nv)); Eh = sp.csr_matrix((np.ones(nh), (np.arange(nh), ih)), shape=(nh, nv))
    # [A  -El^T  Eh^T] [phi    ]   [rhs   ]
    # [El   0     0  ] [lam_lo ] = [lo    ]
    # [Eh   0     0  ] [lam_hi ]   [hi    ]   (symmetric after multiplying the lam_lo column and row by -1: use B = [-El; Eh])
    B = sp.vstack([-El, Eh]).tocsr(); g = np.concatenate([-lo*np.ones(nl), hi*np.ones(nh)]); m = nl + nh
    K = sp.bmat([[A, B.T], [B, None]], format="csc") if m else A.tocsc(); rr = np.concatenate([rhs, g]); it = 0
    if solver == "kkt-lu": s = sla.splu(K).solve(rr)
    elif solver == "kkt-ldl":
        d = 1e-12*abs(A).max(); Kq = sp.bmat([[A, B.T], [B, -d*sp.eye(m)]], format="csc") if m else A.tocsc()
        F = qdldl.Solver(sp.triu(Kq, format="csc")); s = F.solve(rr)
        for _ in range(2): s += F.solve(rr - K @ s)
    elif solver == "kkt-minres":
        DA = A.diagonal(); act = np.concatenate([il, ih]); Pinv = np.concatenate([1/DA, DA[act]]); c = [0]
        s, info = sla.minres(K, rr, x0=np.concatenate([phik, np.zeros(m)]), M=sp.diags(Pinv), rtol=1e-13, maxiter=20000, callback=lambda xk: c.__setitem__(0, c[0]+1)); it = c[0]
    elif solver == "kkt-ppcg":   # projected PCG, constraint preconditioner [diag(A) B^T; B 0] (Gould, Hribar, Nocedal 2001)
        DA = A.diagonal(); act = lo_set | hi_set; x = phik.copy(); x[il] = lo; x[ih] = hi
        Pj = lambda v: np.where(act, 0., v/DA); r = A @ x - rhs; z = Pj(r); p = -z; rz = r @ z; nb = np.sqrt(abs(rz)) + 1e-300   # relative to the initial projected residual
        while it < 20000 and np.sqrt(abs(rz)) > 1e-13*nb and act.sum() < nv:
            Ap = A @ p; al = rz/(p @ Ap); x += al*p; r += al*Ap; z = Pj(r); rz2 = r @ z; it += 1
            if np.sqrt(abs(rz2)) <= 1e-13*nb: break
            p = -z + (rz2/rz)*p; rz = rz2
        lam = A @ x - rhs; s = np.concatenate([x, -lam[il], -lam[ih]])   # B^T mu = rhs - A x -> mu_lo = A x - rhs on il ... sign of the symmetric form
        s[nv:nv+nl] = (A @ x - rhs)[il]; s[nv+nl:] = -(A @ x - rhs)[ih]   # A x - rhs = El^T mu_lo - Eh^T mu_hi
    elif solver == "kkt-uzawa":   # augmented Lagrangian
        rho = 1e2*A.diagonal().max(); Ar = (A + rho*(B.T @ B)).tocsr(); mu = np.zeros(m); x = phik.copy()
        for it in range(1, 200):
            x, _ = sla.cg(Ar, rhs - B.T @ mu + rho*(B.T @ g), x0=x, rtol=1e-13, M=sp.diags(1/Ar.diagonal()))
            v = B @ x - g; mu += rho*v
            if m == 0 or np.abs(v).max() <= 1e-12: break
        s = np.concatenate([x, mu])
    x = s[:nv]; mu = s[nv:]
    lam_lo = np.zeros(nv); lam_hi = np.zeros(nv); lam_lo[il] = mu[:nl]; lam_hi[ih] = mu[nl:]   # A x - rhs = -B^T mu = El^T mu_lo - Eh^T mu_hi
    return x, lam_lo, lam_hi, it
# ---------------- one phase step: primal-dual active set ----------------
def phase_step(phi, eta, rho, lo, hi, solver, stats):
    phin = phi.copy(); phik = phi.copy()
    for kn in range(40):
        t0 = time.perf_counter(); A, r = phase_AR(phik, phik, phin, eta, rho); stats["asm"] += time.perf_counter() - t0
        d = A @ np.ones(nv); rhs = A @ phik - r
        # primal-dual active set with the signed multiplier mu = r (= lambda_lo - lambda_hi at the current iterate), c = d:
        # lower set: mu - d (phi - lo) > 0, upper set: mu - d (phi - hi) < 0 (the solver's test phi - r/d outside [lo, hi])
        tol = 1e-9; stp = phik - r/d
        lo_set = stp <= lo + tol; hi_set = (stp >= hi - tol) & ~lo_set
        stats["act"] = max(stats["act"], int(lo_set.sum() + hi_set.sum()))
        t0 = time.perf_counter()
        x, l_lo, l_hi, it = (solve_elim if solver.startswith("elim") else solve_kkt)(A, rhs, lo_set, hi_set, lo, hi, solver, phik)
        stats["sol"] += time.perf_counter() - t0; stats["kry"] += it; stats["its"] += 1
        x = np.clip(x, lo, hi); dmx = np.abs(x - phik).max(); phik = x
        if dmx < 1e-10:
            A2, r2 = phase_AR(phik, phik, phin, eta, rho); d2 = A2 @ np.ones(nv)
            viol = np.where(phik <= lo + 1e-9, np.maximum(-r2, 0), np.where(phik >= hi - 1e-9, np.maximum(r2, 0), np.abs(r2)))/d2
            stats["kkt"] = max(stats["kkt"], viol.max()); stats["lmin"] = min(stats["lmin"], min(l_lo.min(), l_hi.min()))
            return phik, True
    return phik, False
def run(solver, lo, hi):
    phi = 0.999*np.ones(nv); st = dict(asm=0., sol=0., its=0, kry=0, act=0, kkt=0., lmin=0.); t0 = time.perf_counter(); ok = True
    for s in range(10):
        eta = 1e-3 if s < 5 else 0.1; rho = 10. if s < 5 else 100.
        phi, c = phase_step(phi, eta, rho, lo, hi, solver, st); ok &= c
    st["tot"] = time.perf_counter() - t0
    return phi, st, ok
if __name__ == "__main__":
    ref = np.loadtxt("phi_mode0.txt")
    cases = [("production bounds [0, 1-1e-7]", 0.0, 1 - 1e-7), ("stress bounds [0.3, 0.95]", 0.3, 0.95)]
    solvers = sys.argv[1].split(",") if len(sys.argv) > 1 else ["elim-cg", "elim-lu", "kkt-lu", "kkt-ldl", "kkt-minres", "kkt-ppcg", "kkt-uzawa"]
    if len(sys.argv) > 2: cases = [c for c in cases if c[1] == float(sys.argv[2])]
    names = {"elim-cg": "elimination + CG/Jacobi (current method)", "elim-lu": "elimination + sparse LU", "kkt-lu": "KKT (lam_lo, lam_hi) + sparse LU",
             "kkt-ldl": "KKT + LDL^T (QDLDL)", "kkt-minres": "KKT + MINRES, block-diag. precond.", "kkt-ppcg": "KKT + projected PCG",
             "kkt-uzawa": "KKT + augmented-Lagrangian Uzawa"}
    for tag, lo, hi in cases:
        print(f"\n=== {tag}: 10 phase steps on the production mesh ({nv} nodes)")
        print(f"{'method':42s} {'total s':>8s} {'solve s':>8s} {'assembly s':>10s} {'AS its':>6s} {'Krylov':>7s} {'max active':>10s} {'|phi - first|':>13s} {'KKT viol.':>9s} {'min lambda':>11s}")
        base = None
        for sv in solvers:
            phi, st, ok = run(sv, lo, hi)
            if base is None: base = phi
            print(f"{names[sv]:42s} {st['tot']:8.2f} {st['sol']:8.2f} {st['asm']:10.2f} {st['its']:6d} {st['kry']:7d} {st['act']:10d} {np.abs(phi - base).max():13.1e} {st['kkt']:9.1e} {st['lmin']:11.1e}" + ("" if ok else "  NOT CONVERGED"))
        if lo == 0.0: print(f"Python port against the FreeFEM solver (same steps): max |phi_python - phi_FreeFEM| = {np.abs(base - ref).max():.1e}")
