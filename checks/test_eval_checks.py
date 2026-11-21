#!/usr/bin/env python3
"""Regression test of checks/eval_checks.py: a passing transcript (checks/fixtures: adcheck and loaded-run output and log of a
coarse-mesh self-test, 15 November 2026) and mutations of it that must fail (missing, duplicated, non-finite or malformed records,
failed processes, incomplete runs, failure markers, ERROR lines).  usage: python3 checks/test_eval_checks.py  (exit 1 on a wrong verdict)"""
import os, re, shutil, subprocess, sys, tempfile
here = os.path.dirname(os.path.abspath(__file__)); fx = os.path.join(here, "fixtures")
AD0, RUN0, LOG0 = (open(os.path.join(fx, f)).read() for f in ("ad.out", "t.out", "t_log.csv"))
STEP = 13

def sub_last(text, pat, rep):
    m = list(re.finditer(pat, text, re.M)); assert m, pat
    k = m[-1]; return text[:k.start()] + k.expand(rep) + text[k.end():]
def drop(text, pat, keep=0):
    out, n = [], 0
    for line in text.splitlines(True):
        if re.search(pat, line):
            n += 1
            if n <= keep: out.append(line)
            continue
        out.append(line)
    return "".join(out)

cases = [  # name, expected exit (0 pass, 1 fail), changes
    ("complete transcript", 0, {}),
    ("global tangent record removed", 1, {"run": lambda r: drop(r, r"^TANGENT CHECK \(exact")}),
    ("global tangent error 0.1", 1, {"run": lambda r: re.sub(r"(\(exact Gauss-point tangent\): \|K d - FD\| / \|FD\| = )\S+", r"\g<1>0.1", r)}),
    ("8 of 9 AD records removed", 1, {"ad": lambda t: drop(t, r"^adcheck: output [1-8] ")}),
    ("3 of 4 Mandel records removed", 1, {"run": lambda r: drop(r, r"^FD check Mandel (slip [12]|spin)")}),
    ("NaN in the last AD record", 1, {"ad": lambda t: sub_last(t, r"(adcheck: output 8\s+max \|FD - AD\| )\S+", r"\g<1>nan")}),
    ("-nan in the last AD record", 1, {"ad": lambda t: sub_last(t, r"(adcheck: output 8\s+max \|FD - AD\| )\S+", r"\g<1>-nan")}),
    ("inf denominator in an AD record", 1, {"ad": lambda t: sub_last(t, r"(adcheck: output 4\s+max \|FD - AD\| \S+\s+\(max \|AD\| )[^)]+", r"\g<1>inf")}),
    ("NaN in the last Mandel record", 1, {"run": lambda r: sub_last(r, r"(FD check Mandel spin\s+\(tau = M:B, all elements\): max \|-dpsi/dh - tau\| / max \|tau\| = )\S+", r"\g<1>nan")}),
    ("NaN in a step-scan record", 1, {"run": lambda r: sub_last(r, r"(TANGENT CHECK step scan: FD factor 1e-08\s+\|K d - FD\|/\|FD\| = )\S+", r"\g<1>nan")}),
    ("malformed number in the KWC record", 1, {"run": lambda r: re.sub(r"(KWC force only: \|H d - FD\|/\|FD\| = )\S+", r"\g<1>1.2.3", r)}),
    ("multiplier block B removed", 1, {"run": lambda r: drop(r, r"^TANGENT CHECK multiplier B")}),
    ("global tangent record duplicated", 1, {"run": lambda r: r + [l for l in r.splitlines(True) if l.startswith("TANGENT CHECK (exact")][0]}),
    ("ERROR appended after the diagnostics", 1, {"run": lambda r: r + "ERROR: simulated solver failure after diagnostic output\n"}),
    ("'done.' line missing", 1, {"run": lambda r: drop(r, r"^done\. wall time")}),
    ("loaded run exit status 1", 1, {"run_status": 1}),
    ("adcheck exit status 139", 1, {"ad_status": 139}),
    ("loaded run stopped one step early", 1, {"log": lambda l: "".join(l.splitlines(True)[:-1])}),
    ("NOCONV marker present", 1, {"marker": "_NOCONV.txt"}),
    ("NAN marker present", 1, {"marker": "_NAN.txt"}),
]
bad = 0
for name, want, ch in cases:
    d = tempfile.mkdtemp()
    try:
        ad = ch.get("ad", lambda x: x)(AD0); run = ch.get("run", lambda x: x)(RUN0); log = ch.get("log", lambda x: x)(LOG0)
        if any(k in ch for k in ("ad", "run", "log")):   # a text mutation must change the transcript
            assert (ad, run, log) != (AD0, RUN0, LOG0), "mutation had no effect: " + name
        for f, t in (("ad.out", ad), ("t.out", run), ("t_log.csv", log)): open(os.path.join(d, f), "w").write(t)
        if "marker" in ch: open(os.path.join(d, "t" + ch["marker"]), "w").write("x\n")
        cmd = [sys.executable, os.path.join(here, "eval_checks.py"), "--ad", os.path.join(d, "ad.out"), "--ad-status", str(ch.get("ad_status", 0)),
               "--run", os.path.join(d, "t.out"), "--run-status", str(ch.get("run_status", 0)), "--log", os.path.join(d, "t_log.csv"),
               "--expect-step", str(STEP), "--prefix", os.path.join(d, "t"), "--multiplier", "1"]
        pr = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True); rc = pr.returncode
        why = "; ".join(l.split("  ")[0].strip() for l in pr.stdout.splitlines() if l.rstrip().endswith("FAIL"))
    finally: shutil.rmtree(d)
    ok = (rc == want); bad += (not ok)
    print("%-38s expected %-4s got %-4s %-13s %s" % (name, "pass" if want == 0 else "fail", "pass" if rc == 0 else "fail", "ok" if ok else "WRONG VERDICT", why[:90]))
print("EVALUATOR REGRESSION " + ("PASSED (%d cases)" % len(cases) if bad == 0 else "FAILED (%d wrong verdicts)" % bad)); sys.exit(1 if bad else 0)
