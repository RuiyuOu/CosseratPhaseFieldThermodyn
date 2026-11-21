#!/usr/bin/env python3
# Independence analysis of the matrix check (ASMCHECK every N Newton iterations).
# Input: FEAT lines of a run with -asmcheck 1 (every iteration checked): iteration, step, Newton it, deviation, max T, min phi,
# max |dgamma|, yielding elements, max |dw*|, max |tau|, max |beta|.
import sys, math, numpy as np
from scipy import stats
rows = [l.split()[1:] for l in open(sys.argv[1]) if l.startswith("FEAT ")]
X = np.array([[float(v) for v in r] for r in rows])
K = len(X); N = int(sys.argv[2]) if len(sys.argv) > 2 else 50
it, step, nit, dev, Tm, phm, dg, ny, dw, tau, beta = X.T
print(f"Newton iterations logged: {K} (steps {int(step.min())}..{int(step.max())}); check interval N = {N}")
print(f"matrix check deviation: max {dev.max():.2e}, median {np.median(dev):.2e} (limit 1e-12): no failure in any iteration")

def runs_test(b):
    """Wald-Wolfowitz runs test on a 0/1 sequence: z < 0 means fewer runs than independent (clustering)."""
    b = np.asarray(b, bool); n1 = b.sum(); n0 = len(b) - n1
    if n1 == 0 or n0 == 0: return float("nan"), float("nan")
    r = 1 + np.count_nonzero(b[1:] != b[:-1]); mu = 2*n1*n0/(n1 + n0) + 1
    var = 2*n1*n0*(2*n1*n0 - n1 - n0)/((n1 + n0)**2*(n1 + n0 - 1)); z = (r - mu)/math.sqrt(var)
    return z, 2*stats.norm.sf(abs(z))
def lag1(x):
    x = np.asarray(x, float) - np.mean(x); return float(np.dot(x[1:], x[:-1])/np.dot(x, x)) if np.dot(x, x) > 0 else float("nan")
def episodes(b):
    b = np.asarray(b, bool); L = []; c = 0
    for v in b:
        if v: c += 1
        elif c: L.append(c); c = 0
    if c: L.append(c)
    return L
def miss_sampled(b, N):
    """fraction of the N check phases (checks at iterations p, p+N, ...) that never see a faulty iteration"""
    b = np.asarray(b, bool); return np.mean([not b[p::N].any() for p in range(N)])

print("\n1) The deviation itself (round-off of two assemblies of the same matrix)")
med = np.median(dev); z, p = runs_test(dev > med)
print(f"   lag-1 autocorrelation of log10(dev): {lag1(np.log10(np.maximum(dev, 1e-30))):.2f}; runs test above/below median: z = {z:.2f}, p = {p:.2g}")

print("\n2) State-dependent faults: a branch of the formulas that is wrong only in some states fails exactly when the state enters")
print("   that region. Region = a state variable above a threshold (quantiles chosen so that a fraction f of the iterations fail).")
print(f"   {'state variable':22s} {'f':>6s} {'lag-1':>6s} {'runs z':>7s} {'p':>8s} {'episodes':>8s} {'mean L':>7s} {'max L':>6s} {'P(miss) data':>13s} {'indep.':>9s} {'episode model':>14s}")
feats = {"max T": Tm, "-min phi": -phm, "max |dgamma|": dg, "yielding elements": ny, "max |dw*|": dw}
m = K//N
for name, x in feats.items():
    for f in (0.01, 0.02, 0.05, 0.1, 0.2):
        thr = np.quantile(x, 1 - f); b = x > thr
        if b.sum() == 0 or b.all(): continue
        fe = b.mean(); z, p = runs_test(b); L = episodes(b)
        pind = (1 - fe)**m                                  # independent Bernoulli faults, m checks
        pep = float(np.prod([max(0.0, 1 - l/N) for l in L]))  # each episode of length L missed with probability 1 - L/N
        print(f"   {name:22s} {fe:6.3f} {lag1(b):6.2f} {z:7.2f} {p:8.1e} {len(L):8d} {np.mean(L):7.1f} {max(L):6d} {miss_sampled(b, N):13.3f} {pind:9.2e} {pep:14.3f}")

print("\n3) Faults tied to the position in the Newton iteration (e.g. the first iteration after a start, the first of each step)")
for name, b in {"first iteration of each step": nit == 0, "iteration >= 5 of a step": nit >= 5}.items():
    print(f"   {name:32s} f = {b.mean():.3f}, episodes {len(episodes(b))}, P(miss) with every-{N} sampling = {miss_sampled(b, N):.3f}, independence {(1 - b.mean())**m:.1e}")

print(f"\n4) What independence would allow (m = {m} passed checks of this run; production: about 3000 steps x 7 iterations = 21000, m = 420)")
for mm, Nt in ((m, K), (420, 21000)):
    pu = 1 - 0.05**(1/mm)                                     # 95% upper bound on p, zero failures in mm trials ("rule of three": 3/mm)
    unc = Nt - mm
    print(f"   m = {mm}: 95% upper bound on the per-iteration fault rate {pu:.2e}; expected faulty unchecked iterations at that bound {pu*unc:.1f}; "
          f"P(no fault in the {unc} unchecked | uniform prior) = {(mm + 1)/(mm + unc + 1):.3f}")
