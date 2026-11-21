#!/usr/bin/env python3
"""Pass/fail decision for the self-test (checks/run_tangent_check.sh).

usage: eval_checks.py --ad AD.out --ad-status RC --run RUN.out --run-status RC --log RUN_log.csv --expect-step N
                      --prefix OUTPREFIX [--multiplier 0|1]

A test passes only if all of the following hold (exit status 0, last line SELF-TEST PASSED; otherwise exit 1, SELF-TEST FAILED):
- both FreeFEM processes exited with status 0;
- the loaded run reached the expected final step (last row of RUN_log.csv), printed its 'done.' line exactly once, and left no
  failure marker (<prefix>_NOCONV.txt, _NAN.txt, _TSTOP.txt); neither output contains an ERROR line;
- every expected record is present exactly once: AD outputs 0..8, the plugin-value line, the global tangent, 4 step-scan lines, the
  KWC tangent, the multiplier blocks C, B, E (with --multiplier 1), the Mandel checks slip 0, slip 1, slip 2, spin, the internal force;
- every number used is finite (a malformed or non-finite field fails its test), and each error is within its limit.
Limits are about 100x the values measured on the coarse test mesh (an error of O(1), a wrong derivative, fails clearly). The
finite-increment driving-force comparison is printed for information only (sampling and the finite exponential separate it from
tau = M:B), but its numbers must still be finite.
"""
import argparse, math, os, re, sys

ap = argparse.ArgumentParser()
ap.add_argument("--ad", required=True); ap.add_argument("--ad-status", type=int, required=True)
ap.add_argument("--run", required=True); ap.add_argument("--run-status", type=int, required=True)
ap.add_argument("--log", required=True); ap.add_argument("--expect-step", type=int, required=True)
ap.add_argument("--prefix", required=True); ap.add_argument("--multiplier", type=int, default=1)
a = ap.parse_args()

def read(p):
    try: return open(p, errors="replace").read()
    except OSError: return ""
ad, run = read(a.ad), read(a.run)
rows, fails = [], 0

def result(name, ok, shown, limit="-"):
    global fails
    fails += (not ok); rows.append((name, shown, limit, "pass" if ok else "FAIL"))

def num(tok):
    """float of a printed number; None if malformed or not finite"""
    try: v = float(tok)
    except (TypeError, ValueError): return None
    return v if math.isfinite(v) else None

def once(text, pattern, name):
    """the groups of the single match of pattern; a missing or repeated record fails"""
    m = re.findall(pattern, text, re.M)
    if len(m) != 1:
        result(name, False, "missing" if not m else "%d records" % len(m)); return None
    return m[0] if isinstance(m[0], tuple) else (m[0],)

def limit_test(name, value, limit):
    if value is None: result(name, False, "non-finite", limit); return
    result(name, value <= limit, "%.3g" % value, limit)

T = r"(\S+)"   # any token: the number is validated by num()

# ---- process status and completion ----
result("adcheck process exit status", a.ad_status == 0, str(a.ad_status), "0")
result("loaded-run process exit status", a.run_status == 0, str(a.run_status), "0")
last = None
try:
    for line in open(a.log):
        f0 = line.split(",")[0].strip()
        if f0.isdigit(): last = int(f0)
except OSError: pass
result("loaded run reached the expected final step", last == a.expect_step, str(last), str(a.expect_step))
nd = len(re.findall(r"^done\. wall time", run, re.M))
result("loaded run printed 'done.' once", nd == 1, str(nd), "1")
for mk in ("_NOCONV.txt", "_NAN.txt", "_TSTOP.txt"):
    result("no failure marker %s" % mk, not os.path.exists(a.prefix + mk), "present" if os.path.exists(a.prefix + mk) else "absent", "absent")
ne = len(re.findall(r"ERROR", ad + run))
result("no ERROR line in the outputs", ne == 0, str(ne), "0")

# ---- AD element map: outputs 0..8, each exactly once, each within the limit ----
recs = re.findall(r"^adcheck: output (\S+)\s+max \|FD - AD\| (\S+)\s+\(max \|AD\| (\S+)\)", ad, re.M)
ids = [r[0] for r in recs]; worst, okAD, why = 0.0, True, ""
for i in range(9):
    n = ids.count(str(i))
    if n != 1: okAD, why = False, ("output %d: %d records" % (i, n)); break
if okAD and len(recs) != 9: okAD, why = False, "%d records (expected 9)" % len(recs)
if okAD:
    for _, e, m in recs:
        e, m = num(e), num(m)
        if e is None or m is None or m <= 0: okAD, why = False, "non-finite or zero field"; break
        worst = max(worst, e / m)
if okAD: result("AD element map vs FD, all 9 outputs (worst)", worst <= 1e-5, "%.3g" % worst, 1e-5)
else: result("AD element map vs FD, all 9 outputs (worst)", False, why, 1e-5)
g = once(ad, r"^adcheck: values plugin vs solver max \|diff\| " + T + ";", "plugin vs solver element-map values (absolute)")
if g: limit_test("plugin vs solver element-map values (absolute)", num(g[0]), 1e-9)

# ---- global tangent, step scan, KWC tangent ----
g = once(run, r"^TANGENT CHECK \(exact Gauss-point tangent\): \|K d - FD\| / \|FD\| = " + T, "global tangent vs FD of the residual")
if g: limit_test("global tangent vs FD of the residual", num(g[0]), 1e-5)
sc = re.findall(r"^TANGENT CHECK step scan: FD factor (\S+)\s+\|K d - FD\|/\|FD\| = (\S+)", run, re.M)
fac = [num(f) for f, _ in sc]; val = [num(v) for _, v in sc]
if len(sc) != 4 or len(set(fac)) != 4: result("tangent step scan: 4 FD steps, best", False, "%d records" % len(sc), 1e-5)
elif any(v is None for v in val + fac): result("tangent step scan: 4 FD steps, best", False, "non-finite", 1e-5)
else: result("tangent step scan: 4 FD steps, best", min(val) <= 1e-5, "%.3g" % min(val), 1e-5)
g = once(run, r"^TANGENT CHECK KWC force only: \|H d - FD\|/\|FD\| = " + T, "KWC tangent vs FD")
if g: limit_test("KWC tangent vs FD", num(g[0]), 1e-8)
if a.multiplier:
    for blk in "CBE":
        g = once(run, r"^TANGENT CHECK multiplier " + blk + r" \([^)]*\): \|\S+ \S+ - FD\|/\|FD\| = " + T, "multiplier block %s vs FD" % blk)
        if g: limit_test("multiplier block %s vs FD" % blk, num(g[0]), 1e-5)
elif re.search(r"^TANGENT CHECK multiplier", run, re.M):
    result("no multiplier records without --multiplier 1", False, "present")

# ---- tau = M:B, each generator exactly once ----
for gen in ("slip 0", "slip 1", "slip 2", "spin"):
    g = once(run, r"^FD check Mandel " + gen + r"\s+\(tau = M:B, all elements\): max \|-dpsi/dh - tau\| / max \|tau\| = " + T + r"\s+max \|tau\| = " + T,
             "tau = M:B vs -dpsi/dh, %s" % gen)
    if g:
        e, m = num(g[0]), num(g[1])
        limit_test("tau = M:B vs -dpsi/dh, %s" % gen, e if (m is not None and m > 0) else None, 1e-6)

# ---- internal force ----
g = once(run, r"^FD check internal force: dE/dU\.d = " + T + r"\s+f_int\.d = " + T + r"\s+relative difference = " + T, "internal force vs dE/dU")
if g:
    v = [num(x) for x in g]
    limit_test("internal force vs dE/dU", v[2] if all(x is not None for x in v) else None, 1e-5)

# ---- finite-increment comparison (information; numbers must be finite) ----
fi = re.findall(r"^FD check (slip \d|spin  ) \(finite increment[^)]*\): -dE/dp/area = (\S+)\s+driving force = (\S+)", run, re.M)
if len(fi) != 4: result("finite-increment records (4)", False, "%d records" % len(fi))
for nm, x, y in fi:
    x, y = num(x), num(y)
    if x is None or y is None: result("finite-increment %s numbers finite" % nm.strip(), False, "non-finite")
    else: rows.append(("finite-increment %s (informational)" % nm.strip(), "%.3g" % (abs(x - y) / max(abs(y), 1e-300)), "-", "info"))

w = max(len(r[0]) for r in rows)
for r in rows: print("%-*s  %12s  limit %-6s  %s" % (w, r[0], r[1], r[2], r[3]))
print("SELF-TEST " + ("PASSED" if fails == 0 else "FAILED (%d)" % fails))
sys.exit(1 if fails else 0)
