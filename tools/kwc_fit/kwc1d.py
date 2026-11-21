"""1D KWC tilt boundary with the KWC energy including the quadratic orientation term:
gamma = f0 min int[ (1-phi)^2/2 + a^2 phi'^2/2 + s g(phi) (sqrt(theta'^2+delta^2) - delta) + eps^2/2 h(phi) theta'^2 ] dx,
g(phi) = -2[phi + ln(1-phi)],  h(phi) = 1 (hform='1', default) or phi^2 (hform='phi2'),  theta jumps by dth across the domain.
Alternating minimisation; both sub-problems are convex:
  theta | phi : s g q/sqrt(q^2+delta^2) + eps^2 h q = C in every cell (q = theta' >= 0), C by bisection on int q = dth
                (eps = 0, delta = 0: a jump at the minimum of g)
  phi | theta : convex in u = 1 - phi; L-BFGS-B with the analytic gradient.
Units: lengths um, f0 GPa -> gamma in J/m^2."""
import numpy as np
from scipy.optimize import minimize


def _q_of_C(C, gm, hm, s, eps, delta):
    """per-cell q >= 0 with s gm q/sqrt(q^2+delta^2) + eps^2 hm q = C (monotone in q)"""
    e2h = eps*eps*hm
    if eps == 0.0:
        if delta == 0.0:
            return None
        Cc = np.minimum(C, s*gm*(1 - 1e-15)); return delta*Cc/np.sqrt((s*gm)**2 - Cc*Cc)
    lo = np.zeros_like(gm); hi = np.full_like(gm, C)/np.maximum(e2h, 1e-300)
    for _ in range(80):
        mid = 0.5*(lo + hi)
        F = s*gm*mid/np.sqrt(mid*mid + delta*delta) + e2h*mid
        lo = np.where(F < C, mid, lo); hi = np.where(F < C, hi, mid)
    return 0.5*(lo + hi)


def _theta(gm, hm, h, dth, s, eps, delta):
    if eps == 0.0 and delta == 0.0:
        q = np.zeros_like(gm); i = np.argmin(gm); q[i] = dth/h[i]; return q
    Cmax = s*gm.min()*(1 - 1e-15) if eps == 0.0 else None
    lo, hi = 0.0, (Cmax if Cmax is not None else 1.0)
    if Cmax is None:                                          # grow the bracket until the jump is reached
        while (_q_of_C(hi, gm, hm, s, eps, delta)*h).sum() < dth: hi *= 2
    for _ in range(80):
        C = 0.5*(lo + hi)
        if (_q_of_C(C, gm, hm, s, eps, delta)*h).sum() < dth: lo = C
        else: hi = C
    q = _q_of_C(0.5*(lo + hi), gm, hm, s, eps, delta)
    return q*dth/(q*h).sum()


def gamma(dth, f0, a, s, delta, eps=0.0, L=None, n=1600, iters=80, hform="1"):
    if L is None:
        L = max(8*a, 8*eps, 0.8)
    x = np.linspace(-L, L, n + 1); h = np.diff(x)
    u = 0.7*np.exp(-np.abs(x)/a) + 1e-9

    def E_u(u, q):
        um = 0.5*(u[1:] + u[:-1]); du = np.diff(u); gm = -2*(1 - um + np.log(um)); pm = 1 - um
        r = np.sqrt(q*q + delta*delta) - delta
        hh = pm*pm if hform == "phi2" else 1.0
        e = (h*(0.5*um*um + 0.5*a*a*(du/h)**2 + s*gm*r + 0.5*eps*eps*hh*q*q)).sum()
        dum = h*(um + s*r*(2 - 2/um) - (eps*eps*pm*q*q if hform == "phi2" else 0.0))
        g = np.zeros_like(u); g[:-1] += 0.5*dum; g[1:] += 0.5*dum
        dd = a*a*du/h; g[:-1] -= dd; g[1:] += dd
        return e, g

    Eold = np.inf
    for k in range(iters):
        um = 0.5*(u[1:] + u[:-1]); gm = -2*(1 - um + np.log(um)); hm = (1 - um)**2 if hform == "phi2" else np.ones_like(um)
        q = _theta(gm, hm, h, dth, s, eps, delta)
        r = minimize(E_u, u, args=(q,), jac=True, method="L-BFGS-B", bounds=[(1e-14, 1.0)]*len(u),
                     options=dict(maxiter=5000, ftol=1e-16, gtol=1e-13))
        u = r.x; E = r.fun
        if abs(Eold - E) < 1e-12*abs(E): break
        Eold = E
    c = np.cumsum(q*h)/dth
    w = x[1:][np.searchsorted(c, 0.95)] - x[1:][np.searchsorted(c, 0.05)]
    return f0*E*1e3, dict(phimin=1 - u.max(), theta_width_um=w, outer_its=k + 1)


def gamma_sharp(dth, f0, a, s):
    p = np.linspace(0, 1 - 1e-9, 2000001); g = -2*(p + np.log(1 - p))
    return (f0*(a*(1 - p)**2 + s*dth*g)).min()*1e3
