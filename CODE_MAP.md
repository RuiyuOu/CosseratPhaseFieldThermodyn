# cosserat2d_mob — argument reference and code map (pre-cleanup)

Snapshot of 2026-11-01, updated 2026-11-04 (Section 0). The line numbers in Sections 1–3 refer to these exact files of the snapshot; Section 0 gives the current files, the options added since, and current line numbers in `cosserat2d_mob_stable.edp`:

| file | lines | role |
|---|---|---|
| `cosserat2d_mob.edp` | 1729 | serial solver, the **source of truth** |
| `cosserat2d_mob_mpi.edp` | 1815 | MPI/PETSc solver, **generated** by `python3 mpi/make_mpi.py cosserat2d_mob.edp cosserat2d_mob_mpi.edp` (never edit by hand) |

**Production path** (what the run scripts pass, and what the cleanup keeps):
* `-variant A` (KWC phase field);
* `-tangent consistent`, `-newton full`;
* `-eos vinet`, `-dragrho 3 -vc 2`, `-recov 2 -cd 100`, `-dissip 1 -cdstore 1`;
* `-spinbranch 0.785398 -bckin 1 -predictor vel`;
* `-dir piston` (`ANGLE=confined`) or `-dir shearpiston` (`ANGLE`, default 45);
* linear solver: serial GMRES + ILU(0), umfpack fallback; MPI PETSc GMRES + ASM/ILU, MUMPS fallback.

Status labels used below:
* **KEEP**: on the production path.
* **TOGGLE**: physics or modelling option worth keeping as a switch.
* **DECIDE**: your call before the purge.
* **PURGE**: legacy, superseded or dead (no effect on the production path).
* **DIAG**: diagnostics / tests.

---

## 0. Changes since the snapshot (2026-11-04)

**Files** (the snapshot's `cosserat2d_mob.edp` was split; regenerate every `_mpi` file with `python3 mpi/make_mpi.py <serial> <mpi>` after any edit):

| file | lines | role |
|---|---|---|
| `cosserat2d_mob_stable.edp` | 1909 | stable serial solver (`CODE=stable`, the default), source of truth |
| `cosserat2d_mob_stable_mpi.edp` | 2004 | generated MPI version |
| `cosserat2d_mob_exp.edp` | 2148 | stable + KWC crystal symmetry (`-symfold`, `-symn`, `-kwcfield`, `-kwconeway`, `-kwcp` …; `CODE=exp`, tag `_exp`) |
| `cosserat2d_mob_exp_mpi.edp` | 2248 | generated MPI version |
| `cosserat2d_mob_baek.edp` | 2185 | stable + Baek dual (`-baek 1`, see `BAEK_NOTES.md`); serial only |

**Options added since the snapshot** (line = `cosserat2d_mob_stable.edp`):

| arg | line | solver default | script (default) | meaning | status |
|---|---|---|---|---|---|
| `-gnd` | 728 | `curl` | `GND` (`curl`) | GND density: `curl` = \|G − G_R\|/b with G = Fp Curl Fp (Cermelli–Gurtin); `grad` = \|∇(θ − θ_R)\|/b; `phigrad` = φ²·\|∇(θ − θ_R)\|/b (the old default) | KEEP |
| `-elint`, `-nq` | 125 | `mean`, `3` | `ELINT` (**gauss**), `NQ` (3) | elastic-energy integration: `gauss` = Vinet/Murnaghan finite-strain energy at `-nq` Gauss points per triangle (3 or 6), no quadratic reference energy and no one-point correction (needs `-variant A`, `-eos vinet|murnaghan`, `-mucdual 0`); `mean` = old. New functions `gaussInt`, operators `GD*`/`GX*` (sparse Gauss-point kinematics and their transposes); the tag gets `Q3`/`Q6`. `-nq 3` = the three edge midpoints (FreeFem `qf2pT`). `-gtan` (`GTAN`, default 1): with `-tangent consistent` the Newton matrix is the exact Hessian of the Gauss-point energy (`gaussTangent`: one `int2d(qf2pT)` form with P1dc coefficient fields `gh*`, `gc*`; `gaussTangentSparse` = sparse-product reference, `-gtan 2/3`); `-gtan 0` = old approximate matrix. `-tancheckstep N` runs the `-tancheck` at the first step >= N; `-tancheck 2` (serial) also scans the FD step and checks the KWC force alone. `consistentTangent` now treats a binding spin cap exactly (row of the local system replaced by dq_w = 0) | KEEP |
| `-gndcurl`, `-gndphi2` | 730 | unset | `GNDCURL`, `GNDPHI2` (unset) | legacy; when either is given it overrides `-gnd` (a missing `-gndphi2` then means 1, as before) | PURGE after the old series |
| `-mucdual`, `-dualr` | 46, 48 | `0`, `0` | `MUCDUAL`, `DUALR` | rotational compatibility by multiplier instead of the μ_c penalty | TOGGLE |
| `-kwcediff` | 86 | `1` | `KWCEDIFF` (**0**) | 0 drops the ε²\|∇θ\|² term (tag `_ne`) | KEEP |
| `-kwcreg`, `-kwcgam` | 90 | `0`, `0` (= 1/δ) | `KWCREG`, `KWCGAM` | \|∇θ\| regularisation: 0 = √(p²+δ²) − δ, 1 = Kobayashi–Giga ln cosh | TOGGLE |
| `-etaphi`, `-qphi` | 105 | `0.1`, `100e3` | `ETAPHI` (0.1), `QPHI` (100e3) | φ viscosity at 300 K (GPa ns), Arrhenius energy (J/mol) | KEEP |
| `-etaphimin` | 109 | `0` | `ETAPHIMIN` (**1e-2**) | floor on the Arrhenius factor of η_φ | KEEP |
| `-etastar` | 111 | `1e-3` | `ETASTAR` (**1**) | GB-spin viscosity η*₀ (GPa ns) | KEEP |

**Production defaults changed in the scripts** (from the DRXtest D3 run `pistM3v2RcV_ne_pm1e-2_es1_DRXtest_mpi_a90_up250ms`): `KWCEDIFF` 1 → 0, `ETAPHIMIN` 0 → 1e-2, `ETASTAR` 1e-3 → 1; `RECOV=2`, `CD=100` unchanged. GND: solver and script default `curl` (was φ²∇θ). This is the one solver default that no longer reproduces the old behaviour (user decision); `-gnd phigrad` restores it.

**Tags** mark differences from the *old* script defaults, so new runs never reuse old folders: the default tag is now `pistM3v2GRcVQ3_gc_ne_pm1e-2_es1[_mpi]_a<ANGLE>_up<V>ms` (`Q3`: `-elint gauss`, 2026-11-06). `GND=phigrad KWCEDIFF=1 ETAPHIMIN=0 ETASTAR=1e-3` gives the old `pistM3v2RcV…` tag and physics. The tag logic is the same in `run_piston_mob{,_mpi}.sh`, `run_post_piston_mob{,_mpi}.sbatch` and `resume_all_A.sh` (which previously knew only `G`, `_sm`, `_sf`, `_gc`).

**Wall time (MPI):** `wall=` in the step lines and `-walllimit` now use `mpiWtime()` (elapsed time) instead of rank-0 `clock()` (CPU time, which lagged on contended nodes); the CPU time is printed alongside as `(cpu … s)`.

**Post-processing (2026-11-05):** `post/animate.py` has a new field `rotation` = θ_lattice(t) − θ_lattice(first frame) (diverging scale, ±99.5th percentile of the sampled frames, `ROTLIM=<deg>` fixes it; added to an existing folder without re-rendering the other fields). `post/animate.py` is incremental (per-frame PNGs and `manifest.json`/`limits.json` in `anim_<tag>/`; new frames only, colour limits of the first call kept; `FORCE=1` redoes all, `RELIM=1` also new limits). `post/grain_post.py` runs `post/grain_count.py` on every frame (relaxed frame + `<tag>.pvd`) into `grains_<tag>/grain_history.csv` (also incremental; grain parameters `GRAD MERGE MINAREA SYM LO MID HI` in `state.json`, a change re-analyses all) and writes `grain_history.png`, the animations `grain_theta`, `grain_territory`, `grain_misorientation` (PNGs + GIF in `anim_<tag>/`, fixed colour scales, territories coloured by mean theta; third argument of `grain_post.py`) plus `grain_table_all.csv` (one row per grain and frame) and `boundary_table_all.csv` (one row per neighbour pair and frame) and full tables and three figures (`_theta.png`, `_territory.png`, `_misorientation.png`) of the first and last frame. Both are called by `run_post_piston_mob{,_mpi}.sbatch`; `package.sh` gets the grain folder as 5th argument and copies it to `deliver_<tag>/grains/`.

**Phase-field Arrhenius factor:** `macro phArr()` = max(exp(Q_φ/R (1/T − 1/T0)), r_min) is used in the φ solve and in the dissipation; `make_mpi.py` patches the `etaPhi*phArr` form.

---

## 1. Command-line arguments (serial and MPI)

"Script" is the environment variable in `run_piston_mob.sh` / `run_piston_mob_mpi.sh` and its default there; `—` means not exposed by the scripts.
Solver defaults reproduce the old behaviour (exception since 2026-11-04: `-gnd curl`); the production values live in the scripts.

### 1.1 Run control, loading, output

| arg | line | solver default | script (default) | meaning | status |
|---|---|---|---|---|---|
| `-variant` | 29 | `F` | 6th positional (`A`) | `A` = KWC phase field + GB spin; `F` = Forest plastic curvature | DECIDE (F: purge?) |
| `-j` | 30 | `1e-2` | fixed `0.01` | micro-inertia j (µm²) | KEEP |
| `-amp` | 31 | `0.06` | fixed `1` | traction amplitude (GPa), traction modes only | DECIDE (with traction modes) |
| `-dir` | 32 | `shear` | `ANGLE` → `piston` / `shearpiston` | loading: `shear`, `comp`, `combo` (tractions), `piston` (confined piston), `shearpiston` (angled plate) | KEEP piston, shearpiston; DECIDE traction modes |
| `-shear` | 33 | `0` | — | `combo`: extra shear traction (GPa) | DECIDE (traction) |
| `-tend` | 34 | `15` | 2nd positional (`30`) | end time (ns) | KEEP |
| `-pistonup` | 35 | `0.25` | V·cos(ANGLE) | normal plate / piston speed (µm/ns) | KEEP |
| `-pistonut` | 36 | `0` | V·sin(ANGLE) | tangential plate speed (µm/ns), shearpiston only | KEEP |
| `-ptr` | 37 | `0.5` | `RAMP` (3rd positional, 0.5) | ramp time (ns) | KEEP |
| `-pth` | 38 | `2.0` | `HOLD` (4th positional, 2.0) | hold time (ns) | KEEP |
| `-dt` | 39 | `0.01` | `DT` (0.01) | time step (ns) | KEEP |
| `-muc` | 40 | `75.4` | fixed | Cosserat coupling modulus (GPa) | KEEP |
| `-out` | 41 | `run_<variant>` | `out_<tag>/<tag>` | output prefix | KEEP |
| `-nrelax` | 42 | `100` | — | A: zero-load relaxation steps | KEEP |
| `-vtk` | 43 | `1` | — | write VTU frames | KEEP |
| `-vtkevery` | 44 | `0` | 5th positional (10, scaled by 0.01/DT) | extra frames every N steps | KEEP |
| `-maxit` | 45 | `40` | `MAXIT` (40) | global Newton iterations | KEEP |
| `-verbose` | 46 | `0` | — | 1: per-iteration and timing lines; 2: per-element local iterations | DIAG |
| `-tol` | 47 | `1e-7` | — | Newton increment tolerance du | KEEP |
| `-stepsmax` | 48 | `1e8` | — | stop after this step | KEEP |
| `-tscale` | 49 | `1` | — | traction envelope time scale | DECIDE (traction) |
| `-restart` | 50 | `""` | checkpoint of the tag | resume from a checkpoint | KEEP |
| `-tstop` | 51 | `1e30` | `TSTOP` (1400 K) | graceful stop above max T | KEEP |
| `-walllimit` | 52 | `1e30` | — | checkpoint and stop after N s of wall time | KEEP |
| `-ckevery` | 53 | `100` | 50·(0.01/DT) | regular checkpoint interval | KEEP |
| `-ckroll` | 111 | `0` | `CKROLL` (20) | rolling pre-crash checkpoints `_ckroll{0,1}` | KEEP |
| `-relaxfix` | 1446 | `1` | `RELAXFIX` (1) | apply the end-of-relaxation reset to a relaxed checkpoint | KEEP (drop the 0 branch) |
| `-trace` | 1139 | `1` | — | `_stage.txt` breadcrumbs | KEEP |
| `-diag` | 665 | `0` | `DIAG` (1) | per-step worst-element `_diag.csv` | KEEP |

### 1.2 Material / KWC phase field (variant A)

| arg | line | solver default | script | meaning | status |
|---|---|---|---|---|---|
| `-f0` | 73 | `0.0115` | fixed `0.005944` | KWC energy scale (GPa) | KEEP (set script value as default) |
| `-kwca` | 73 | `0.062` | fixed `0.10` | gradient coefficient a (µm), sets the core width | KEEP |
| `-kwcs` | 73 | `0.150` | fixed `1.05` | coupling s | KEEP |
| `-kwce` | 74 | `0.200` | fixed `0.004` | ε (θ-gradient) | KEEP |
| `-kwcd` | 74 | `0.005` | fixed `0.005` | regularisation δ of \|∇θ\| | KEEP |
| `-etastar` | 77 | `1e-3` | `ETASTAR` (1 since 2026-11-04) | GB spin friction η*₀ (GPa ns) | KEEP |
| `-cd` | 77 | `10` | `CD` (100 with RECOV=2) | sweeping-recovery constant C_D | KEEP |
| `-recov` | 666 | `0` | `RECOV` (2) | recovery: 0 = C_D\|dφ\|, 1 = ReLU, 2 = Ask 2019 ReLU·tanh | KEEP 2; DECIDE 0/1 |
| `-ca` | 667 | `3.16227766` | `CA` | C_A (µm) for recov 2 | KEEP |
| `-gndphi2` | 664 | `1` | `GNDPHI2` (1) | φ² weight on the solver GND. **Superseded by `-gnd` (Section 0)** | legacy |

### 1.3 Elastic law

| arg | line | solver default | script | meaning | status |
|---|---|---|---|---|---|
| `-eos` | 82 | `svk` | `EOS` (`vinet`) | volumetric law: `svk` (old), `murnaghan`, `vinet` | TOGGLE (you asked to keep murnaghan/svk) |
| `-k0` | 84 | `0` (= cubic K = 137.1) | `K0` | bulk modulus (GPa) | KEEP |
| `-kp` | 85 | `5.0` | `KP` | K′ | KEEP |
| `-voltan` | 86 | `1` | `VOLTAN` (1) | (K(J) − K₀)⁺ correction in the **approx** tangent | PURGE (only acts in the approx tangent) |

### 1.4 Slip law, drag, mobile density

| arg | line | solver default | script | meaning | status |
|---|---|---|---|---|---|
| `-dragrho` | 98 | `0` | `DRAGRHO` (3) | ρ_m in the drag: 0 = ρ_S+ρ_G, 1 = ρ_S, 2 = Ma–Roters, 3 = evolved ρ_M | KEEP 3; DECIDE 0–2 |
| `-vc` | 116 | `0` | `VC` (2 km/s) | relativistic drag limit | KEEP |
| `-rhomsat` | 119 | `100` | `RHOMSAT` (100) | ρ_M saturation (µm⁻²) | KEEP (see the rate-ceiling note in README) |
| `-c2m` | 119 | `10` | `C2M` | ρ_M evolution rate | KEEP |
| `-rhofcut` | 119 | `1e4` | `RHOFCUT` | forest cutoff | KEEP |
| `-rhommin` | 119 | `1` | `RHOMMIN` | ρ_M floor | KEEP |
| `-c1mr`, `-c2mr`, `-c3mr` | 118 | `0.1`, `2.0`, `1.0` | — | Ma–Roters constants (dragrho 2 only) | PURGE with dragrho 2 |
| `-softmin` | 115 | `0.01` | `SOFTMIN` | floor of the thermal softening factor | KEEP |

### 1.5 Robustness and thermodynamics

| arg | line | solver default | script | meaning | status |
|---|---|---|---|---|---|
| `-spinbranch` | 107 | `0` | `SPINBRANCH` (π/4) | principal-root bound \|dw*\| ≤ b (spin update + polish) | KEEP (make π/4 the solver default) |
| `-spinmax` | 101 | `0` | `SPINMAX` (0) | physical bound on \|dw*\| per step | PURGE (superseded by spinbranch) |
| `-onfail` | 109 | `retry` | `ONFAIL` (retry) | unconverged step: retry / stop / commit | KEEP retry/stop, PURGE commit |
| `-maxitretry` | 110 | `80` | `MAXITRETRY` | retry iterations | KEEP |
| `-dissip` | 113 | `0` | `DISSIP` (1) | heat source τ_{n+1}Δγ (1) vs midpoint (0) | KEEP 1; DECIDE whether `_thermo.csv` keeps the midpoint audit |
| `-cdstore` | 114 | `0` | `CDSTORE` (1) | storage cap Z dρ ≤ cdstore·Σ\|τΔγ\| | KEEP |
| `-bckin` | 917 | `0` | `BCKIN` (1) | exact kinematics on the prescribed dofs | KEEP 1, PURGE 0 (0 is the bug) |
| `-predictor` | 918 | `avg` | `PREDICTOR` (vel) | Newton start: avg / vel / none | KEEP vel, PURGE avg/none |

### 1.6 Tangent, Newton, linear algebra

| arg | line | solver default | script | meaning | status |
|---|---|---|---|---|---|
| `-tangent` | 269 | `approx` | `TANGENT` (consistent) | consistent element tangent vs plastic condensation | KEEP consistent, PURGE approx |
| `-fdk`, `-fdq` | 271 | `1e-6`, `1e-8` | `FDK`, `FDQ` | FD steps of the element map in k and q | KEEP |
| `-polish` | 272 | `20` | `POLISH` | coupled local Newton iterations | KEEP |
| `-nlocal` | 646 | `30` | `NLOCAL` | robust local slip/spin iterations | KEEP |
| `-newton` | 1323 | `modified` | `NEWTON` (full) | full Newton + line search vs Aitken fixed point | KEEP full, PURGE modified |
| `-tolres` | 1325 | `1e-6` | `TOLRES` | relative residual tolerance | KEEP |
| `-resfloor` | 1325 | `1e-8` | `RESFLOOR` | residual floor | KEEP |
| `-lsmax` | 1326 | `6` | `LSMAX` | line-search halvings | KEEP |
| `-nref` | 647 | `10` | `NREF` (20) | tangent refresh interval (modified Newton only) | PURGE (also drop it from the tag) |
| `-tanevery` | 306 | `1` | `TANEVERY` (1) | full Newton: rebuild every n iterations | PURGE (no gain measured) |
| `-linsolver` | 300 | `umfpack` | `LINSOLVER` (gmres) | serial: GMRES + ILU or umfpack; ignored by MPI | KEEP gmres (umfpack as fallback only) |
| `-lintol` | 302 | `1e-10` | `LINTOL` | linear tolerance (also the PETSc rtol) | KEEP |
| `-krylov`, `-linmaxit` | 303 | `100`, `1000` | — | serial GMRES restart and max iterations | KEEP |
| `-precevery` | 308 | `0` | `PRECEVERY` (0) | ILU rebuild interval | PURGE (always rebuilt with the tangent in full Newton) |
| `-kcvarf` | 304 | `1` | — | consistent tangent by one varf (1) or sparse products (0) | KEEP 1, PURGE 0 |
| `-slippred` | 305 | `0` | `SLIPPRED` (0) | start the local solve from the previous step's increments | PURGE (no gain measured) |
| `-rfr` | 1509 | `4` | — | tangent refresh every n relaxation steps (approx tangent) | PURGE (full Newton rebuilds anyway) |
| `-tancheck` | 1345 | `0` | — | one FD check of the assembled tangent | DIAG (keep) |
| `-check` | 1699 | `0` | — | energy FD checks of driving forces and internal force (serial only) | DIAG (DECIDE) |

### 1.7 MPI-only (cosserat2d_mob_mpi.edp)

| arg | MPI line | default | script | meaning |
|---|---|---|---|---|
| `-npexpect` | 34 | `0` | `NP` | stop with LAUNCHER ERROR unless mpisize = N |
| `-petsc` | 347 | `-ksp_type gmres -ksp_gmres_restart 200 -ksp_max_it 2000 -ksp_rtol <lintol> -pc_type asm -pc_asm_overlap 1 -sub_pc_type ilu` | `PETSC` | (u, θ) Newton system |
| `-petscfall` | 348 | `-ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps` | `PETSCFALL` | redo of a solve that misses lintol |
| `-petscp1` | 349 | `-ksp_type gmres -ksp_rtol 1e-13 -ksp_max_it 2000 -pc_type asm -sub_pc_type lu` | `PETSCP1` | phase-field and heat systems |

### 1.8 Script-only variables (not solver args)

| variable | default | meaning |
|---|---|---|
| `ANGLE` | 45 | plate angle from −x (0 compression … 90 shear); `confined` = old piston. Tag `_a<ANGLE>` |
| `RUNTAG` | — | free label in the tag |
| `GND` | curl | solver `-gnd` (curl / grad / phigrad); tag `G…_gc` / `G` / none |
| `NP`, `MEM`, `MPIRUN`, `SHARETAG` | 8, 24G, auto, 0 | MPI ranks / memory / launcher / reuse the serial tag |
| `N` (run_all, resume_all) | 10 (MPI) / 20 | chained 24 h array links |
| `MPI`, `TBLOW`, `FORCE` | 0, 5000, 0 | resume_all_A.sh: series selector, blow-up threshold, override the running check |

---

## 2. Serial solver `cosserat2d_mob.edp`: block map

| lines | block | status |
|---|---|---|
| 1–19 | header comment | KEEP (rewrite at cleanup) |
| 20–21 | `load "iovtk"`, `load "IncompleteCholesky"` | KEEP |
| 23–28 | `arg()` helper | KEEP |
| 29–53 | run-control arguments (§1.1) | KEEP / DECIDE traction args |
| 56–77 | material constants, KWC parameters, η*, C_D | KEEP |
| 79–86 | elastic-law options (eos, k0, kp, voltan) | KEEP / PURGE voltan |
| 87–98 | mobile-drag documentation, `-dragrho` | KEEP |
| 99–111 | robustness options (spinmax, spinbranch, onfail, maxitretry, ckroll) | KEEP / PURGE spinmax |
| 112–119 | thermodynamics (dissip, cdstore, softmin), vc, Ma–Roters, ρ_M constants | KEEP / PURGE c*mr |
| 120–123 | Schmid factors, slip directions | KEEP |
| 125–135 | mesh, grain angles, fespaces Wh = [P2,P2,P1], Vh1, Ph; T, φ, θ0s, θR | KEEP |
| 136–196 | rotated cubic stiffness, Schmid tensors, G and K matrices per grain → per-element arrays (K11..K23 used by F only) | KEEP (F parts: DECIDE) |
| 197–198 | P0 stiffness fields c11..c66 | KEEP |
| 200–217 | history arrays (slip, ρ_S, ρ_G, ρ_M, w*, F^p, trial F^p, s⊗n) | KEEP (kp1/kp2/dva..dvc: F only) |
| 219–225 | A: smoothing of the initial orientation field θ0s | KEEP |
| 227–236 | operator macros (symv, rotv, sig, dbl), curvK | KEEP |
| 237–242 | `vK` linear elastic + Cosserat (+ F curvature) operator | KEEP (Ktrue in the residual) |
| 243–248 | `vKa` approximate tangent | PURGE (0 refs) |
| 249–250 | mass `vMu`, micro-inertia `vMt` | KEEP |
| 251–260 | BC labels LBa..LBf (piston, shearpiston), `vBC` | KEEP |
| 261–264 | Mu, Mt, Mall, Ktrue, Aeff, Kpl | KEEP (Kpl: PURGE) |
| 265–267 | `usePlTan`, `vR`, `vDs` | PURGE (vR, vDs: 0 refs) |
| 268–274 | `-tangent`, fdk/fdq/polish, P0 coefficient fields cD00..cD44 | KEEP |
| 275–277 | `vKc` consistent tangent alone | PURGE (0 refs) |
| 278–282 | `vKvC` KWC + consistent | PURGE (0 refs) |
| 283–293 | **`vAall`**: whole matrix in one assembly (elastic + Cosserat + inertia + KWC tangent + consistent element tangent + Dirichlet) | KEEP (production matrix) |
| 294 | `Kcons` | PURGE (kcvarf 0 path) |
| 295–311 | linear-algebra options; Kconst, ILU matrices, `PreconILU` | KEEP gmres/ILU; PURGE precevery/slippred/tanevery |
| 312–314 | `vKv` KWC + EOS volumetric tangent (approx path) | PURGE |
| 315–336 | `buildA`: Kconst (first call), vAall (consistent) or Kconst + Kv (+ Kcons / −Kpl), GMRES/ILU set-up | KEEP the vAall + ILU branch; PURGE Kconst/Kv/Kcons/Kpl branch |
| 337–349 | loads: Tload, traction direction, `envelope` | DECIDE (traction) |
| 350–354 | `vF` traction + plastic RHS varf | PURGE (0 refs) |
| 355–359 | KWC macros gtx/gty/gphi/mKfac | KEEP |
| 360–362 | `vKWC`, `vBCr` | PURGE (0 refs) |
| 363–386 | precomputed RHS operators: B1..B3, Bw (plastic distortion), Bk1/Bk2 (F), Fun (unit traction), Bg1/Bg2 (KWC couple stress), BC11..BCth (finite-strain correction) | KEEP (Bk1/Bk2: F; Fun: traction) |
| 387 | `rpl` | KEEP |
| 389–393 | heat: `vTh`, ATh (umfpack), Dsrc, `vTr` | KEEP |
| 395–405 | local-update work arrays; GB-spin elimination arrays | KEEP |
| 406–416 | `effDrag`: relativistic drag η_eff and dη/d\|τ\| | KEEP |
| 417–435 | `flow`: rate law \|γ̇\|(τ) and derivative | KEEP |
| 437–438 | `dragBound`, `sgn` | KEEP |
| 440–463 | `solve3`: 3×3 per-element Cramer solve (robust slip Newton) | KEEP |
| 465–470 | `clampv` | KEEP |
| 472–485 | finite-deformation work arrays | KEEP |
| 486–513 | `trialFp`: F^p_trial = exp(dL^p) F^p_n (closed-form exponential; the series branch and `-exactexp` were removed) | KEEP |
| 514–560 | `eosStress`: vol/dev split, Vinet / Murnaghan, cubic SVK on the isochoric strain | KEEP (TOGGLE eos) |
| 561–623 | `finiteStress`: Ū, Green strain, S, Π, P_F, ψ; quadratic counterparts; correction dP; Mandel stress, τ* | KEEP |
| 625–633 | `taus`: resolved shear stresses τ_a, τ* | KEEP |
| 635–644 | `nus`: F plastic-curvature driving forces | DECIDE (F) |
| 646–647 | `-nlocal`, `-nref` | KEEP nlocal; PURGE nref |
| 648–657 | per-element backtracking | KEEP |
| 658–662 | `setHp`: plastic fields for the assembly | KEEP |
| 663–668 | locres; `-gndphi2`, `-diag`, `-recov`, `-ca`; centroids xcen/ycen | KEEP |
| 669–675 | `kin`: element-mean kinematics | KEEP |
| 676–685 | consistent-tangent work arrays (DSF, DSL, DTQ, JJ, BB, DCT) | KEEP |
| 686–704 | `solve4`: in-place elimination of nt 4×4 systems | KEEP |
| 706–717 | `addIn`, `outI`: FD perturbation helpers | KEEP |
| 719–736 | `fdMap`: FD of the element map s(k,q), τ(k,q) | KEEP |
| 738–746 | `buildJqq`: exact local Jacobian J_qq | KEEP |
| 748–756 | `localRes`: local residuals | KEEP |
| 758–790 | `polishLocal`: coupled slip + spin Newton with backtracking and the spinbranch clamp | KEEP |
| 791–860 | `localUpdate`: robust slip loop (+ spin step, drag bound), F curvature loop, polish, Hp | KEEP (F loop: DECIDE) |
| 863–883 | `strength`: τ_c(ρ, T), ρ_m (dragrho modes), η_d, rate ceiling, spin friction | KEEP (dragrho 0–2 branches: DECIDE) |
| 885–893 | state vectors U, V, Acc; `plasticRHS` | KEEP |
| 894–902 | BC masks: bcm/mask, piston (pistm, pistk), shear plate (pistmT), rollers (confm, confk) | KEEP |
| 903–910 | `pistonDisp` | KEEP |
| 911–918 | bcm/mask (all Dirichlet dofs), `-bckin`, `-predictor` | KEEP |
| 919–949 | `pistonVel`, `pistonAcc`, `pistonDispU/VelU/AccU` (plate profile) | KEEP (merge into one unit-profile set at cleanup) |
| 950–955 | `fixBCkin`: exact kinematics on prescribed dofs | KEEP |
| 956–970 | KWC element arrays, areas, `updatePhi2`, `kwcStress` | KEEP |
| 971–980 | `residual` | KEEP |
| 981–988 | Rk (slip-strain operators) for the approx tangent | PURGE |
| 989–997 | Kk (element-mean kinematics operators) | KEEP? (only the kcvarf 0 path uses Kk/KkT → PURGE with it) |
| 1000–1057 | `consistentTangent`: FD map, J_qq, dq/dk, D_tot − D_lin → cD fields (1016–1043); sparse-product path 1044–1055 | KEEP; PURGE 1044–1055 (kcvarf 0) |
| 1058–1074 | cp*, Iv*, N*, W*, Dm* arrays for the plastic condensation | PURGE |
| 1075–1134 | `plasticTangent` (approx tangent) | PURGE |
| 1135–1145 | `-trace`, `mark` | KEEP |
| 1147–1173 | `nanE`, `nanV` | KEEP |
| 1174–1199 | `refreshTangent`: KWC tangent H, consistent tangent (or plastic condensation + voltan), buildA | KEEP the consistent branch; PURGE the else branch |
| 1200–1222 | energies: elastic, stored, kinetic (u, θ), heat | KEEP |
| 1224–1233 | output files, headers (log, diag, thermo) | KEEP |
| 1234–1246 | P0→P1 helper (F only), cos/sin θ0, gid, pvd, outbase | KEEP (toP1: F) |
| 1247–1298 | `writeVTK` (7 P1 + 31 P0 fields) | KEEP |
| 1299–1317 | time-loop set-up: initial strength/local/tangent, step counts, work arrays | KEEP (omg/rprev/dr: PURGE with modified Newton) |
| 1318–1329 | `-newton`, `-tolres`, `-resfloor`, `-lsmax` | KEEP |
| 1330–1343 | `linSolve` (GMRES with residual check, umfpack fallback) | KEEP |
| 1345–1359 | `-tancheck`, `checkTangent` | DIAG |
| 1360–1413 | `newtonSolve`: local update → tangent → residual → convergence/divergence tests → linear solve → line search (full) / Aitken (modified) | KEEP; PURGE the Aitken branch |
| 1414–1415 | banner | KEEP |
| 1416–1425 | `saveCkF`, `saveCk` (+ `.rhoM` side file) | KEEP |
| 1426–1455 | restart: read checkpoint, ρ_M side file, pvd, relaxfix reset, bckin reset, initial strength/local/tangent | KEEP |
| 1456–1493 | `_meta.json` | KEEP |
| 1494 | F: t0 frame | DECIDE (F) |
| 1495–1697 | **time loop** (rows below) | |
| 1496–1510 | step start: relax flag, load, inertial RHS rin, stage, φ² weights, strength (1509: `-rfr` refresh) | KEEP (PURGE 1509) |
| 1511–1520 | predictor (avg / vel / none) + exact piston / plate displacement | KEEP vel |
| 1521–1548 | global Newton; retry with a fresh tangent; clean NOCONV stop; timing line | KEEP |
| 1549–1550 | final local update; scheduled refresh (modified Newton only) | KEEP 1549, PURGE 1550 |
| 1551–1558 | Newmark update, fixBCkin, external work | KEEP |
| 1559–1575 | dissipation (midpoint and implicit), GB spin, F curvature dissipation | KEEP (F part: DECIDE) |
| 1576–1577 | KMT storage + dynamic recovery | KEEP |
| 1578–1604 | phase-field Newton (backward Euler, log coupling), sweeping recovery (recov 0/1/2), φ dissipation | KEEP |
| 1605–1612 | Clausius–Duhem storage cap; D bookkeeping | KEEP |
| 1613–1618 | commit slip, w*, ρ_S, F^p, slip paths; NaN guards | KEEP (slip predictor line: PURGE) |
| 1619–1624 | ρ_M evolution (dragrho 3) | KEEP |
| 1625–1634 | GND: F from plastic curvature, A from ∇(θ − θ_R) · φ² (now: `-gnd` curl / grad / phigrad, Section 0) | KEEP (F: DECIDE) |
| 1635–1640 | heat equation | KEEP |
| 1641–1650 | end of relaxation: relaxed checkpoint and reset | KEEP |
| 1651–1682 | per-step output: energies, log, thermo, diag, console | KEEP |
| 1683–1685 | VTU frames (fixed times, every vtkevery) | KEEP |
| 1686–1696 | NaN during commit, checkpoints, rolling checkpoints, TSTOP, wall limit | KEEP |
| 1698 | end-of-run summary | KEEP |
| 1699–1729 | `-check` energy FD checks | DIAG (DECIDE) |

## 3. MPI solver `cosserat2d_mob_mpi.edp`: MPI-specific blocks

Everything not listed is the serial code shifted by the inserted lines (+86 overall); the patch script keys every change on exact text, not line numbers.

| MPI lines | block |
|---|---|
| 21–23 | `load "PETSc"`, `macro dimension`, `include "macro_ddm.idp"` (replaces IncompleteCholesky) |
| 33–39 | `-npexpect` launcher guard |
| 136–141 | global mesh ThG, `buildDmesh(Th)` → overlapping subdomain, n2o, nt/ntG |
| 147–173 | global spaces WhG/Vh1G/PhG; Mats AeffM (vector), Ap1 (φ), AThM (heat), Ap0; ownership ownE/ownW/own1; interpolation RW/R1; reductions gsum/gmax/gmin/gdotW/gsumE; gathers gathE/gathW/gath1; PETSc option prefixer `pfx` |
| 260–264 | θ0s smoothing solved on the global mesh, restricted |
| 345–349 | linMode = 0; `-petsc`, `-petscfall`, `-petscp1` |
| 361–378 | `buildA` (local matrix → `AeffM`, `set` with the `el_` prefix) |
| 431–432 | heat Mat AThM (`ht_` prefix) |
| 799–860 | local loops with global termination tests (gmax in localUpdate, polishLocal) |
| 1188–1215 | `nanE`/`nanV` with collective detection; per-rank NaN reports |
| 1247–1266 | energies with owned-element weights + gsum |
| 1271–1276 | log / diag / thermo opened on rank 0 only (/dev/null elsewhere) |
| 1331–1348 | `writeVTK`: gather of 9 P1 + 31 P0 fields, savevtk on ThG by rank 0 |
| 1387–1399 | `linSolve`: PETSc solve, owned-dof residual check, MUMPS fallback |
| 1401–1415 | `checkTangent` with a global direction (same as serial) |
| 1416–1470 | `newtonSolve` with global norms (gdotW, gmax) |
| 1472–1487 | `saveCkF`: gather to the global numbering, rank-0 write in the serial format |
| 1491–1500 | restart: global checkpoint read on every rank, restricted with RW/R1/n2o |
| 1656–1666 | phase field: varf + Ap1 PETSc solve (`ph_`), global dmx test |
| 1714 | heat solve `AThM^-1` |
| 1750–1765 | diag: gathered worst elements on rank 0 |
| 1785 | `-check` disabled (serial only) |
| — | every `cout` guarded by `if (mpirank == 0)` (regex in make_mpi.py) |

## 4. Proposed purge (on the production path these have no effect)

**Dead code (0 references), remove outright:** `vKa` 243–248, `vR`/`vDs` 266–267, `vKc` 275–277, `vKvC` 278–282, `vF` 350–354, `vKWC` 360–361, `vBCr` 362.

**Approximate tangent and modified Newton (superseded by consistent + full):**
* the tangent itself: `plasticTangent` + its arrays 1058–1134; `Kpl`, `usePlTan`; `Rk*` 981–988;
* the assembly path: `vKv` 312–314; the Kconst/Kv/Kcons branch of `buildA`; the else branch of `refreshTangent`; `voltan`;
* modified Newton: the Aitken branch of `newtonSolve` (omg, rprev, dr); `nref` and its tag part; the scheduled refresh 1550; `-rfr` 1509;
* the kcvarf 0 sparse path 1044–1055 with `Kk*` 989–997 and `Kcons`;
* the `-tangent`, `-newton`, `-kcvarf` switches themselves.

**Tuning flags without measured benefit:** `-tanevery`, `-precevery`, `-slippred` (+ dpaP..dwstP), `-krylov` (keep if the serial GMRES stays).

**Superseded or buggy options:** `-spinmax`, `-bckin 0`, `-predictor avg/none`, `-onfail commit`, `-relaxfix 0`.

**Decide first:**
* Variant F (curvature plasticity: `nus`, dva..dvc, kq, Bk1/Bk2, K11..K23, toP1, F branches in localUpdate, dissipation and GND). The consistent tangent exists for A only.
* Traction load modes (`shear`/`comp`/`combo`, `-amp`, `-shear`, `-tscale`, `envelope`, Fun, Wext): the plate modes use none of them.
* `-dragrho` 0/1/2 and the Ma–Roters constants.
* `-recov` 0/1.
* `-eos svk`/`murnaghan` (you asked to keep them as toggles).
* `-dissip 0` (the midpoint audit in `_thermo.csv`).
* The `-check` block.

**After the purge:** make the production values the solver defaults, so the scripts pass only physics and run control. Rename the `pistonVel/Acc/DispU` set into one unit-profile function. Regenerate the MPI file with make_mpi.py; its patches key on exact text, so any patched line that moves fails loudly and must be updated together.

## Reviewer-driven safeguards (2026-11-07), `-safe 1`
Options block just before `elintName` in `cosserat2d_mob_stable.edp`: `-stagres`, `-jmin` (`GV1..GV4`, `jBound(U)`: Bernstein lower bound of det F, called at the start of `residual()`), `-lsfinal` and stagnation handling in `newtonSolve`, `-phistrict` (phase-field loop in the commit; stops with `_NOCONV.txt`), `-heatfloor` (SSD floor inside the stored-energy increment), `-mobend` (`mobUpdate()`), `-platework` (`Rfull`, `cdof`, `reactionNow()`, `Wplate`, side file `.plate`, CSV columns `W_plate,audit_full`), `-gndcg` (in `curlFpG`), atomic checkpoints (`saveCkA`). Serial wrappers `gmaxS/gminS/gsumS/gdotS/gsumEs/isRoot` are identities and are replaced by `mpi/make_mpi.py`. See `docs/cosserat_kwc_stable.md` section 13 for the disposition of each reviewer item.

Second review (same date): `saveCkA` writes immutable generations `<ck>.g<step>_<k>` and publishes `<ck>.manifest` (restart resolves it right after `-restart` is read); `curlFpG` has separate `-gndcg` and legacy branches; the phase loop has the active set (`penU`, `penL`, `vphR`, `vphD`, `syncPen`) and the KKT test (`-phikkt`); `plateLog` owns the plate-work ledger and baseline (`resR`, `plateBaseline`). `-exactexp` and the series evaluation of the plastic exponential were removed (`trialFp` is the closed form only).

Third review: the phase step is a primal-dual active-set (Lagrange-multiplier) solve with exact elimination of the bound-active nodes (`vphA`, `vphB`, `vphL`, `pmask`, `xpin`; no penalty); `.plate` is the tagged `PLATEv2` format with a legacy reader; `-expreject` (in `residual`, flag `qRawMaxLast` from `trialFp`) rejects plastic increments with |q| > 50; the launchers and `recover_ck.sh` follow `<ck>.manifest`.

PETSc notes for the multiplier phase step: `docs/cosserat_kwc_stable.md`, section 15 (implementation notes): elimination on the local matrices with synchronised flags; active nodes set exactly to their bound after the (iterative) solve, bound tolerance 1e-9 in the KKT test; `set(Am, solver = sparsesolver)` needed in the serial code.
Plate sidecar `.plate`: three layouts (`PLATEv2` with `END`: 7+N tokens; patch-2 untagged: 5+N; patch-1 untagged: 3+N) validated on the token stream (array-length token = N and exact total count) before any vector read; anything else takes the explicit unavailable/new-baseline path.

`recover_ck.sh` is a checked transaction (unique `.g<step>_rb_<time>_<pid>` generation, validation, one `mv -f` of the manifest, markers cleared only after commit, exit 1 on failure); `resume_all_A.sh` skips a run whose recovery failed. `ckstep` in the launchers validates manifest, generation and step before any log trimming.

**tau / temperature / mismatch statistics (2026-11-08):** `post/tau_stats.py RUN_DIR [...] --window T0 T1 --mark 6 --log <tag>_log.csv` reads every `<tag>_s<n>.vtu` of a run (times from the `.pvd`; triangle cells only, the boundary line cells are skipped) and writes a per-snapshot CSV (all / bulk / core selections by cell-mean phi), `tau_summary.csv` and `tau_vs_time.png` (four panels: max|tau| statistics and tau_c; per-system tau; T and pressure; Cosserat mismatch beta). Quantities: `tau_1..3` (Mandel resolved shear, GPa, signed; statistics of |tau|), `tau_c`, T (cell mean), pressure -(sig_11+sig_22+sig_33)/3, the largest rise of mean/max T (heating jump), and beta = phi_e - dtheta with phi_e = atan2(G21-G12, G11+G22), G = (I+grad u) Fp^-1 (diagnostic of the penalty residual; `max_abs_dtheta_deg` of the solver log is the rotation itself). Warns if the series ends before `--tend` (30 ns) or has gaps > 0.5 ns. Called by `run_post_piston_mob{,_mpi}.sbatch` into `tau_<tag>/`; `package.sh` copies it to `deliver_<tag>/tau/` (6th argument).
