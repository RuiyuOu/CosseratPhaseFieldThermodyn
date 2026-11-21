#!/usr/bin/env python3
"""apply_prod80.py -- production freeze 8.0 of CosseratPhaseFieldThermodyn (in place, run in the repository root).

  python3 apply_prod80.py            dry run: checks the tree and prints every action, changes nothing
  python3 apply_prod80.py --apply    does it

What it does (nothing under out_*/, deliver_*/, anim_*/, remesh*/, diag_*/ is touched):
  * refuses unless the solver sources are byte-identical to the 7.5 tree that was reviewed (sha256 below);
  * tag.sh: default -kwctan 3 (the kt3pcg set-up); the PHYS string is otherwise unchanged (checked for several environments);
  * run.sh: one line per start in <tag>.out saying whether the solver files match SOLVER_SHA256;
  * VERSION 8.0, SOLVER_SHA256, checks/verify_solver.sh, tools/README.md, checks/phase_solver/README.md, post/CHANGELOG.md;
  * docs: README.md, CODE_MAP.md, cosserat_kwc_production.md (Sections 6.12, 6.13 from r8_fastasm/DOCS_r8.md), v2/README.md;
  * moves: fitting tools -> tools/kwc_fit/, debug generator -> tools/debug/, r8 benchmark scripts -> checks/phase_solver/,
    tests/gnd_cg_test.edp -> checks/;
  * removes the applied patch scripts (apply_*.sh/.py, r8_fastasm/), the duplicate root compare_fields.py, and the *.pre_* backups
    left by earlier patch scripts.
Every file it changes, moves or removes is first saved in _preprod80_backup_<date>.tgz (restore: tar xzf <that file>).
Each text edit must match its anchor exactly the stated number of times, or the script stops before changing anything.
"""
import hashlib, os, re, subprocess, sys, tarfile, time, glob, shutil

APPLY = "--apply" in sys.argv
ROOT = os.getcwd()
def die(m): sys.exit("STOPPED: " + m + "\n(nothing was changed)" if not CHANGED else "STOPPED: " + m)
CHANGED = False
def sha(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()
def rd(p): return open(p, encoding="utf-8").read()

# ---------------------------------------------------------------- 0. the tree must be the reviewed 7.5 tree
for f in ["VERSION", "cosserat_kwc.edp", "tag.sh", "run.sh", "README.md"]:
    if not os.path.isfile(f): die(f"run this in the repository root ({f} missing)")
ver = rd("VERSION").strip()
if ver == "8.0": sys.exit("VERSION is already 8.0: nothing to do")
if ver != "7.5": die(f"VERSION is {ver}, expected 7.5")
FROZEN = {   # solver sources and generated files of the reviewed tree; none of them is changed by this script
    "cosserat_kwc.edp":        "da262e0fb3bf41487a7a9b5fd62dac339d423489edc9483eb8b141f4eeac14f3",
    "v2/make_v2.py":           "6590d2f7d796d7b2e1cd0262d1f22fce76d5dac589982b83f6cf843b1ab11007",
    "mpi/make_mpi.py":         "dececff2c15d55da90888ea4d350a60e00c8ede78b36351b72c386fa4b1ccdd8",
    "plugin/kwc2d.cpp":        "1c1870d79fc6e4d87a76207975aec9aeef25b53a8ebf656e1530cb5e71d2f956",
    "plugin/dual.hpp":         "7dbd3fac13d8b7dbc23883c8860e8af5c392e3e75e9696a5a97364e47f5d8370",
    "cosserat_kwc_v2.edp":     "5a870a5ac6f2fc4791b9bdc29fcb8a497c1d788b0981a0571a9e55d6c5791aa7",
    "cosserat_kwc_v2_mpi.edp": "ad6a7704a9f511672dec6c51f3b1e21ad94fc95bef4181ea835ddb31ac5dbda8",
}
for f, h in FROZEN.items():
    if not os.path.isfile(f): die(f"{f} missing")
    if sha(f) != h: die(f"{f} differs from the reviewed 7.5 tree (sha256 {sha(f)[:12]}...): send it before freezing")
EDITED_IN = {"tag.sh": "e2c6a60c187493e1dfa0a1306b19183b5ec4d4e0cedc1f796cd9b3da4006ae57",
             "run.sh": "a52cc1da22fb105e719624f28adc2982b8ef7c33e27559e103e3b1e9377f426c"}
for f, h in EDITED_IN.items():
    if sha(f) != h: die(f"{f} differs from the reviewed 7.5 tree: send it before freezing")
DOCS_R8 = "r8_fastasm/DOCS_r8.md"
if not os.path.isfile(DOCS_R8) or sha(DOCS_R8) != "3888e7f282c4b9817fef301530cd0f9e17ec8181c96302af0a618e704a52e34b":
    die(f"{DOCS_R8} missing or changed (its Sections 6.12/6.13 go into the production notes)")

# ---------------------------------------------------------------- 1. text edits, prepared in memory
NEW = {}   # path -> new text
def edit(path, old, new, count=1):
    s = NEW.get(path, rd(path) if os.path.exists(path) else None)
    if s is None: die(f"{path} missing")
    n = s.count(old)
    if n != count: die(f"{path}: anchor found {n} times, expected {count}:\n  {old[:160]!r}")
    NEW[path] = s.replace(old, new)

# tag.sh ---------------------------------------------------------------
edit("tag.sh", """#      Both set-ups run the same solver (cosserat_kwc_v2*.edp) with the same numerics and parameters (Cosserat constraint by a multiplier
#      per element, mu_c 754 its augmentation; delta, tolerances, safeguards), except the micro-inertia j and therefore the theta integrator:""",
"""#      Production (version 8.0, frozen): the defaults below are the kt3pcg set-up (exact primal-dual KWC tangent -kwctan 3, Cosserat
#      multiplier -almc 4; the projected-CG phase solve and the fast assembly are solver defaults outside this key), so a run without
#      physics variables is the production run.
#      Both set-ups run the same solver (cosserat_kwc_v2*.edp) with the same numerics and parameters (Cosserat constraint by a multiplier
#      per element, mu_c 754 its augmentation; delta, tolerances, safeguards), except the micro-inertia j and therefore the theta integrator:""")
# the fall-backs after J, MUC and KWCD were never used (all three are always set above): write the variables, same string
edit("tag.sh", 'PHYS="-nrelax ${NRELAX:-100} -j ${J:-0.01} -muc ${MUC:-75.4} -f0 ${F0:-0.005944} -kwca ${KWCA:-0.10} -kwcs ${KWCS:-1.05} -kwcd ${KWCD:-0.005}"',
               'PHYS="-nrelax ${NRELAX:-100} -j $J -muc $MUC -f0 ${F0:-0.005944} -kwca ${KWCA:-0.10} -kwcs ${KWCS:-1.05} -kwcd $KWCD"')
edit("tag.sh", 'PHYS="$PHYS -kwctan ${KWCTAN:-2} -almc $ALMC"   # KWC tangent: 2 = primal-dual (symmetrised), 3 = exact primal-dual Newton; Cosserat coupling: 4 = multiplier, 0 = penalty',
               'PHYS="$PHYS -kwctan ${KWCTAN:-3} -almc $ALMC"   # KWC tangent: 3 = exact primal-dual Newton (production), 2 = primal-dual symmetrised; Cosserat coupling: 4 = multiplier, 0 = penalty')

# run.sh ---------------------------------------------------------------
edit("run.sh", """echo "[$(date +%T)] code version $(cat VERSION 2>/dev/null || echo unknown)" >> $D/$B.out
""", """echo "[$(date +%T)] code version $(cat VERSION 2>/dev/null || echo unknown)" >> $D/$B.out
# frozen solver (version 8.0): the solver sources, the generated files and tag.sh against SOLVER_SHA256 (checks/verify_solver.sh)
if [ -s SOLVER_SHA256 ]; then
  if sha256sum -c --quiet SOLVER_SHA256 >/dev/null 2>&1; then echo "[$(date +%T)] solver frozen: files match SOLVER_SHA256" >> $D/$B.out
  else echo "[$(date +%T)] solver MODIFIED: $(sha256sum -c --quiet SOLVER_SHA256 2>/dev/null | tr '\\n' ' ')" >> $D/$B.out; fi
fi
""")

# README.md ------------------------------------------------------------
edit("README.md", "| 7.5 | 17 Nov 2026: GND measure", """| 7.6 | 18-19 Nov 2026: fast assembly of every staggered field by the plugin on fixed patterns (phase matrix and residual, element means of φ² and g(φ), phase dissipation, GND curl, heat right-hand side), compared with the FreeFEM forms on the first, every 50th and each first-in-regime call (`-fastcheck`, `-fastchecktol`); `FASTASM=0` restores the forms. Phase-field linear solve by projected CG (CG + Jacobi on the eliminated system; `-phsolver cg`, MPI `-petscph`) with a direct fallback. Storage-cap counter in the safeguard line and `<out>_storecap.csv`. Fix: a reduction evaluated on rank 0 only hung MPI runs at the end of the relaxation. |
| 8.0 | 20 Nov 2026: **production freeze.** The defaults are the kt3pcg set-up (`-kwctan 3`, exact primal–dual KWC tangent; projected-CG phase solve; fast assembly): a run without physics variables is the production run. Solver files checksummed (`SOLVER_SHA256`, `bash checks/verify_solver.sh`; every `.out` says whether they match). Applied patch scripts removed; fitting and debug tools moved to `tools/`; post-processing versioned on its own (`post/CHANGELOG.md`). |
| 7.5 | 17 Nov 2026: GND measure""")
edit("README.md", "the KWC term uses the primal–dual (Chan–Golub–Mulet) Newton linearisation;",
     "the KWC term uses the exact reduced Newton tangent of the primal–dual (Chan–Golub–Mulet) system, unsymmetrised (`-kwctan 3`);")
edit("README.md", "  - tangent matrix assembled directly by the plugin on a fixed sparsity pattern;",
     "  - the tangent matrix and every matrix and vector of the staggered fields assembled directly by the plugin on fixed sparsity\n    patterns, compared with the FreeFEM forms every 50th call;")
edit("README.md", "  - staggered phase-field solve (primal-dual active set) and heat solve.",
     "  - staggered phase-field solve: primal–dual active set, each iteration solved by projected CG (CG + Jacobi on the free nodes;\n    Gould, Hribar and Nocedal 2001); heat solve.")
edit("README.md", "\n## Installation\n", """
## Production freeze (version 8.0)

The solver is frozen at the kt3pcg set-up. Frozen files (sha256 in `SOLVER_SHA256`): `cosserat_kwc.edp`, `v2/make_v2.py`,
`mpi/make_mpi.py`, `plugin/kwc2d.cpp`, `plugin/dual.hpp`, the generated `cosserat_kwc_v2.edp` and `cosserat_kwc_v2_mpi.edp`, and `tag.sh`
(the parameters). `bash checks/verify_solver.sh` checks them, and every start of `run.sh` writes `solver frozen` or `solver MODIFIED` into
the run's `.out` file.

- **Production run:** no physics variables. `ANGLE=90 ARRAY=3-5 bash submit_speeds.sh` gives the tags `a90_up150ms_v2` … `a90_up400ms_v2`.
  The kt3pcg folders (`RUNTAG=kt3pcg KWCTAN=3`) have the same physics under their own tag; keep using `RUNTAG=kt3pcg` for them.
- **What still changes:** post-processing (`post/`, its own log in `post/CHANGELOG.md`), the documentation, and the submit and helper
  scripts. None of these can change a solution.
- **Not production:** `remesh/` and `remeshtest/` are the mesh-convergence studies; they hold their own solver copies from earlier
  versions. `tools/` holds the KWC-parameter fit and the debug-build generator (`tools/README.md`).

## Installation
""")
edit("README.md", "| KWC | `-f0` 0.005944, `-kwca` 0.10, `-kwcs` 1.05, `-kwcd` 0.005 |",
     "| KWC | `-f0` 0.005944, `-kwca` 0.10, `-kwcs` 1.05, `-kwcd` 5e-4 |")
edit("README.md", "| Cosserat, elasticity | `-muc` 75.4, `-j` 0.01, `-kp` 5.0, `-k0` 0, `-softmin` 0.01 |",
     "| Cosserat, elasticity | `-muc` 754 (augmentation of the multiplier), `-j` 0 (`OLD=1`: 0.01), `-kp` 5.0, `-k0` 0, `-softmin` 0.01 |\n| tangent, coupling | `-kwctan` 3 (exact primal–dual), `-almc` 4 (multiplier) |")
edit("README.md", "| output | `-diag` 1, `-ckroll` 20 |\n",
     "| output | `-diag` 1, `-ckroll` 20 |\n\nSolver options outside the settings key (recorded per start in `<tag>.xargs`): fast assembly `-asm 1 -fastasm 1` (`FASTASM=0` for the\nFreeFEM forms), `-fastcheck` 50, `-fastchecktol` 1e-8, `-phsolver` cg, exact slip law, `LOCALTOL` 1e-10, `BETAMAX` 0.785398; MPI solver\noptions `PETSC`, `PETSCFALL`, `PETSCP1`, `PETSCPH` (phase field, default CG + Jacobi, rtol 1e-13).\n")
edit("README.md", "`checks/probe_mumps.edp` test a FreeFem module's MPI/PETSc and MUMPS.\n",
     "`checks/probe_mumps.edp` test a FreeFem module's MPI/PETSc and MUMPS. `bash checks/verify_solver.sh` checks the frozen solver files\nagainst `SOLVER_SHA256`. `checks/phase_solver/` holds the benchmark scripts behind the phase-solver comparison (production notes, 6.13).\n")
edit("README.md", "| `checks/` | self-tests and run comparison |",
     "| `checks/` | self-tests, run comparison, `verify_solver.sh`, `phase_solver/` benchmarks |\n| `SOLVER_SHA256`, `VERSION` | checksums of the frozen solver files; the code version |\n| `tools/` | KWC-parameter fit to Tschopp et al. (`kwc_fit/`), debug-build generator (`debug/`); not used by the runs |\n| `remesh/`, `remeshtest/` | mesh-convergence studies (own solver copies from earlier versions; not production) |")

# CODE_MAP.md ----------------------------------------------------------
edit("CODE_MAP.md", "Where things are in the code. How to run it: `README.md`. What each output column and field means: `DATA_FIELDS.md`.\n",
     "Where things are in the code. How to run it: `README.md`. What each output column and field means: `DATA_FIELDS.md`.\n\n"
     "Version 8.0 is the production freeze: the files listed in `SOLVER_SHA256` are not edited any more (`checks/verify_solver.sh`).\n")
edit("CODE_MAP.md", "production runs (both set-ups) pass `-kwctan 2 -almc 4 -asm 1`",
     "production runs (both set-ups) pass `-kwctan 3 -almc 4 -asm 1`")
edit("CODE_MAP.md", "\n## MPI generator `mpi/make_mpi.py`\n", """
## Fast assembly of the staggered fields and the phase solve (version 7.6, in `cosserat_kwc.edp`)

- options `-fastasm` (`fastAsm`), `-fastcheck`, `-fastchecktol`; `fastDue` (which calls are compared), `phRegime` (first call with the
  Arrhenius floor / a node on a bound), `relRow` (row-wise relative comparison), `absMat`, `fastStop` (`<out>_FASTCHECK.txt`);
- phase field: plugin `k2dPhSetup` (pattern), `k2dPhAsm` (matrix and residual, 7-point rule), `k2dPhElem` (element means of φ², g(φ),
  phase dissipation); diagonal $A\\mathbf 1$ and right-hand side $A\\phi_k-R$ by algebra;
- GND: `CurlA`, `CurlB` (constant P0 → P1 operators), `curlFast`; heat right-hand side by constant operators;
- phase linear solve: `-phsolver cg` (`phCG`), `phPCG` (CG + Jacobi on the eliminated system, direct fallback, `nPhCGfail`), a failed
  phase step repeated once with the previous solver (`nPhFall`); MPI: PETSc CG + Jacobi, option `-petscph` (`mpi/make_mpi.py`);
- storage cap: `nStoreCapEl`, `nStoreCapSteps`, `EstoreCapCum`, `<out>_storecap.csv`.

## MPI generator `mpi/make_mpi.py`
""")
edit("CODE_MAP.md", "\n## Checks `checks/`\n", """
## Tools `tools/` (not used by the runs)

- `kwc_fit/`: `kwc_fit_tschopp.py`, `tschopp_grab.py`, `kwc1d.py`, `README_tschopp.md` (fit of the KWC constants to Tschopp et al. 2015).
- `debug/`: `make_dbg.py` (debug build of the MPI solver, `cosserat_kwc_v2_dbg_mpi.edp`); `diag_step.sh` (repository root) reruns
  steps after a checkpoint with per-iteration output.

## Checks `checks/`
""")
edit("CODE_MAP.md", "- `probe_mpi.edp`, `probe_mumps.edp`: installation probes for a FreeFem module.\n",
     "- `probe_mpi.edp`, `probe_mumps.edp`: installation probes for a FreeFem module.\n"
     "- `verify_solver.sh`: `sha256sum -c SOLVER_SHA256` (the frozen solver files).\n"
     "- `phase_solver/`: `kkt2.py` (Python port of the phase step, KKT solver comparison), `bench.edp` (FreeFEM timing of the linear solves),\n"
     "  `ind.py` (independence analysis of the periodic matrix check); `gnd_cg_test.edp`: test of the Cermelli–Gurtin GND measure.\n")

# cosserat_kwc_production.md -------------------------------------------
P = "cosserat_kwc_production.md"
edit(P, "- both production set-ups (`-kwctan 2`): the **primal–dual** linearisation",
        "- `-kwctan 2`: the **primal–dual** linearisation, symmetrised")
edit(P, "**Exact primal–dual Newton tangent (`-kwctan 3`, opt in).**",
        "**Exact primal–dual Newton tangent (`-kwctan 3`, production since version 8.0).**")
edit(P, "not the tangent. Production keeps `-kwctan 2`.",
        "not the tangent. Production uses `-kwctan 3` since version 8.0 (the kt3pcg runs); `-kwctan 2` remains available.")
edit(P, "the free nodes solve the reduced system.",
        "the free nodes solve the reduced system, by projected CG (Section 6.13).")
d = rd(DOCS_R8); i = d.index("## 6.12 Fast assembly")
sec = d[i:].replace("(`patch_fastasm_r8.py`, version 7.2)", "(version 7.6)") \
          .replace("(serial: `LinearCG`;", "(serial: the loop `phPCG`, stopping at $\\|b-Ax\\|_\\infty\\le10^{-13}\\|b\\|_\\infty$;") \
          .replace("Tested\nserially (coarse mesh, and with CG capped at one iteration to exercise the fallback); the MPI path is generated, not yet run.",
                   "Tested\nserially (coarse mesh, and with CG capped at one iteration to exercise the fallback) and in the 16-rank kt3pcg production runs.")
if "not yet run" in sec or "LinearCG" in sec: die("DOCS_r8.md: the MPI-status sentence was not replaced")
edit(P, "\n# 7. Safeguards\n", "\n" + sec.rstrip() + "\n\n# 7. Safeguards\n")
edit(P, "| one rank per physical core | `run.sh` (`HWGUARD=0` to override) | 10 |",
        "| fast assembly compared with the forms (first, every N-th and first-in-regime call), stop above the tolerance | `-fastcheck` 50, `-fastchecktol` $10^{-8}$ | 6.12 |\n"
        "| phase CG: unconverged solve redone by the direct solver; failed phase step repeated once with the previous solver | built in (counted) | 6.13 |\n"
        "| storage cap: elements and steps where the stored energy is capped by the dissipation (`<out>_storecap.csv`) | built in (counted) | 5.5 |\n"
        "| one rank per physical core | `run.sh` (`HWGUARD=0` to override) | 10 |")
edit(P, "| `-petscp1` | GMRES + ASM/LU | $\\phi$ and $T$ |",
        "| `-petscph` | CG + Jacobi, rtol $10^{-13}$ | $\\phi$ (projected CG, Section 6.13) |\n| `-petscp1` | GMRES + ASM/LU | $T$; $\\phi$ when a phase step is repeated with the previous solver |")
edit(P, "  -diag 1 -ckroll 20 -kwctan 2 -almc 4", "  -diag 1 -ckroll 20 -kwctan 3 -almc 4")
edit(P, "| Numerics | primal–dual KWC tangent; direct assembly; exact slip law;",
        "| Numerics | exact primal–dual KWC tangent (unsymmetrised); direct assembly of all fields; projected-CG phase solve; exact slip law;")

# v2/README.md ---------------------------------------------------------
edit("v2/README.md", "| `-kwctan 0\\|2\\|3` | 2 (`KWCTAN`) |", "| `-kwctan 0\\|2\\|3` | 3 (`KWCTAN`, since 8.0) |")

# new files ------------------------------------------------------------
NEW["VERSION"] = "8.0\n"
NEW["checks/verify_solver.sh"] = """#!/bin/bash
# verify_solver.sh -- the frozen solver files (version 8.0) against SOLVER_SHA256; exit 0 if all match.
cd "$(dirname "${BASH_SOURCE[0]}")/.." || exit 1
sha256sum -c SOLVER_SHA256 && echo "solver frozen: all files match SOLVER_SHA256"
"""
NEW["tools/README.md"] = """# tools/ (not used by the runs)

- `kwc_fit/`: fit of the KWC constants to the grain-boundary energies of Tschopp et al. (2015): `kwc_fit_tschopp.py`, `tschopp_grab.py`
  (digitised data), `kwc1d.py` (1D KWC profiles), `README_tschopp.md`. Run from this folder; outputs go to the current directory.
- `debug/`: `make_dbg.py` writes a debug build of the MPI solver (`python3 tools/debug/make_dbg.py cosserat_kwc_v2_mpi.edp
  tools/debug/cosserat_kwc_v2_dbg_mpi.edp`, from the repository root). For rerunning steps of a run use `diag_step.sh` (root).
  `make_dbg.py` was written for version 7.3; some of its anchors no longer match the 8.0 solver (it stops and names the one).
"""
NEW["checks/phase_solver/README.md"] = """# Phase-field solver benchmarks (version 7.6)

Records of the comparison in the production notes, Section 6.13. Run from the repository root (they read `mesh.msh` and
`grain_angles_rad.txt`).
- `kkt2.py`: Python port of the phase step (matches the solver to 1.3e-11) with six solvers of the KKT system with explicit multipliers.
- `bench.edp`: FreeFEM timing of the linear solves of the active-set step (needs the standalone phase plugin it loads; mode in the last
  argument).
- `ind.py`: independence analysis of the periodic matrix check from the FEAT lines of a run with `-asmcheck 1`.
"""
NEW["post/CHANGELOG.md"] = """# Post-processing changelog

Post-processing is versioned separately from the solver, which is frozen at version 8.0 (`SOLVER_SHA256`). Nothing in `post/` or in
`post.sbatch`, `post_all.sh` and `grain_summary.sbatch` can change a solution; record every change of a figure, a statistic or a default
here, with the date, so that a figure in the manuscript can be traced to the version that drew it.

| Date | Change |
|---|---|
| 2026-11-20 | Baseline at the solver freeze: `animate.py`, `grain_post.py` (grain recovery, misorientation; `MINAREA`), `grain_count.py`, `grain_summary.py`, `tau_stats.py`, `wavepost.py`, `sliplaw_compare.py`, `compare_fields.py`, `package.sh`. |
"""

# ---------------------------------------------------------------- 2. moves and removals
MOVES = [("kwc1d.py", "tools/kwc_fit/kwc1d.py"), ("kwc_fit_tschopp.py", "tools/kwc_fit/kwc_fit_tschopp.py"),
         ("tschopp_grab.py", "tools/kwc_fit/tschopp_grab.py"), ("README_tschopp.md", "tools/kwc_fit/README_tschopp.md"),
         ("make_dbg.py", "tools/debug/make_dbg.py"), ("cosserat_kwc_v2_dbg_mpi.edp", "tools/debug/cosserat_kwc_v2_dbg_mpi.edp"),
         ("r8_fastasm/tests/kkt2.py", "checks/phase_solver/kkt2.py"), ("r8_fastasm/tests/bench.edp", "checks/phase_solver/bench.edp"),
         ("r8_fastasm/tests/ind.py", "checks/phase_solver/ind.py"), ("tests/gnd_cg_test.edp", "checks/gnd_cg_test.edp")]
REMOVE = ["apply_lattice_tests.sh", "apply_storecap.py", "apply_v2.sh", "compare_fields.py", "r8_fastasm", "tests"]
for a, b in MOVES:
    if not os.path.exists(a): die(f"{a} missing (expected in the 7.5 tree)")
    if os.path.exists(b): die(f"{b} exists already")
if sha("compare_fields.py") != sha("post/compare_fields.py"): die("compare_fields.py differs from post/compare_fields.py: keep it or merge first")
for r in REMOVE:
    if not os.path.exists(r): die(f"{r} missing (expected in the 7.5 tree)")
# leftover files of r8_fastasm/ and tests/ that are not moved must be known ones
known = {"r8_fastasm/" + f for f in ["APPLY.md", "apply_r8.py", "apply_r8c.py", "apply_r8d.py", "DOCS_r8.md", "patch_fastasm_r8.py",
                                     "r8_plugin_block.cpp", "tests/kkt2.py", "tests/bench.edp", "tests/ind.py", "tests/feat_full.txt"]} | {"tests/gnd_cg_test.edp"}
extra = [os.path.join(dp, f) for top in ["r8_fastasm", "tests"] for dp, _, fs in os.walk(top) for f in fs]
extra = [p for p in extra if p not in known and not p.endswith(".pyc")]
if extra: die("unexpected files in r8_fastasm/ or tests/ (move them first): " + ", ".join(extra))
SKIP = re.compile(r"^(out_|deliver_|anim_|remesh|diag_|grain_summary|logs|\.git)")
PRE = sorted(p for d in [".", "v2", "mpi", "plugin", "post", "checks", "cluster", "mesh"] for p in glob.glob(os.path.join(d, "*.pre_*")) if os.path.isfile(p))
PRE = [os.path.normpath(p) for p in PRE]

# ---------------------------------------------------------------- 3. report / dry run
print(f"tree OK: version 7.5, solver files identical to the reviewed tree ({len(FROZEN)} files)")
print("edit:  " + ", ".join(sorted(p for p in NEW if os.path.exists(p))))
print("new:   " + ", ".join(sorted(p for p in NEW if not os.path.exists(p))) + ", SOLVER_SHA256")
for a, b in MOVES: print(f"move:  {a} -> {b}")
print("remove: " + ", ".join(REMOVE) + (f", {len(PRE)} backup files *.pre_*: " + " ".join(PRE) if PRE else ""))
top = sorted(f for f in os.listdir(".") if os.path.isfile(f))
ODD = [f for f in top if re.search(r"\.(zip|tgz|tar|gz|orig|bak)$", f)]
if ODD: print("left alone (decide yourself): " + " ".join(ODD))
if not APPLY:
    print("\ndry run: nothing changed. Run again with --apply."); sys.exit(0)

# ---------------------------------------------------------------- 4. apply
OLD_TAG = rd("tag.sh")
bk = f"_preprod80_backup_{time.strftime('%Y%m%d_%H%M%S')}.tgz"
with tarfile.open(bk, "w:gz") as t:
    for p in sorted(set([p for p in NEW if os.path.exists(p)] + [a for a, _ in MOVES] + REMOVE + PRE)): t.add(p)
print(f"backup: {bk}")
CHANGED = True
for p, s in NEW.items():
    os.makedirs(os.path.dirname(p) or ".", exist_ok=True)
    st = os.stat(p).st_mode if os.path.exists(p) else None
    with open(p + ".tmp80", "w", encoding="utf-8") as f: f.write(s)
    os.replace(p + ".tmp80", p)
    if st is not None: os.chmod(p, st)
os.chmod("checks/verify_solver.sh", 0o755)
for a, b in MOVES:
    os.makedirs(os.path.dirname(b), exist_ok=True); shutil.move(a, b)
for r in REMOVE:
    shutil.rmtree(r) if os.path.isdir(r) else os.remove(r)
for p in PRE: os.remove(p)
FILES = list(FROZEN) + ["tag.sh"]
with open("SOLVER_SHA256", "w") as f:
    for p in FILES: f.write(f"{sha(p)}  {p}\n")

# ---------------------------------------------------------------- 5. checks after the change
EXPECT_TAG = "53f45306a0ecdd5b37a8b10be35907f61f2a326c1958be7aeba9709d7b221686"
if EXPECT_TAG[0] != "@" and sha("tag.sh") != EXPECT_TAG: print("WARNING: tag.sh differs from the tested result")
def phys(tagtext, env):
    open(".tag80.sh", "w").write(tagtext)
    e = {k: v for k, v in os.environ.items() if k not in ("OLD", "INERTIA", "KWCTAN", "J", "MUC", "KWCD", "RUNTAG", "ALMC")}
    e.update(env); e["UPMS"] = "150"
    out = subprocess.run(["bash", "-c", "source ./.tag80.sh >/dev/null && echo \"$PHYS\""], env=e, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True).stdout.strip()
    os.remove(".tag80.sh"); return out
bad = 0
for env in [{}, {"OLD": "1"}, {"KWCTAN": "3"}, {"KWCTAN": "2"}, {"MUC": "100"}, {"J": "0.02"}, {"KWCD": "1e-3"}, {"OLD": "1", "KWCTAN": "3"}]:
    a, b = phys(OLD_TAG, env), phys(rd("tag.sh"), env)
    if "KWCTAN" not in env: a = a.replace("-kwctan 2", "-kwctan 3")
    if a != b or not a: bad += 1; print(f"PHYS DIFFERS for {env}:\n  old {a}\n  new {b}")
print("PHYS check: " + ("identical to 7.5 apart from the -kwctan default (8 environments)" if not bad else f"{bad} differences (see above)"))
prod = phys(rd("tag.sh"), {})
for a in sorted(glob.glob("out_*_v2_kt3pcg/*.args")):
    k = re.sub(r" -dt [^ ]+ -pistonup .*$", "", re.sub(r" -fd[kq] [^ ]+", "", rd(a).strip()))
    print(f"{a}: " + ("same physics as the 8.0 defaults" if k == prod else "DIFFERS from the 8.0 defaults:\n  run  " + k + "\n  8.0  " + prod))
r = subprocess.run(["sha256sum", "-c", "--quiet", "SOLVER_SHA256"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
print("SOLVER_SHA256: " + ("all files match" if r.returncode == 0 else "MISMATCH " + r.stdout))
print("done: version 8.0. Undo: tar xzf " + bk + " && rm -rf tools checks/phase_solver checks/gnd_cg_test.edp checks/verify_solver.sh SOLVER_SHA256 post/CHANGELOG.md")
