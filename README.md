# HPC bundle: dynamic Cosserat-KWC 15 ns, variants A and F
Contents: cosserat2d.edp (finite-deformation solver), mesh_A.msh, mesh_F.msh, run_dyn.sh, run_dyn_{A,F}.sbatch, run_post.sbatch,
submit_all.sh, post/ (wavepost.py, animate.py, package.sh, SMOOTHING.md).
Use: unzip, cd into the folder, `bash submit_all.sh [N]` submits a chain of N (default 4) restart-safe jobs per variant: each resumes from the last checkpoint, so no manual resubmission is needed. Manual: `sbatch run_dyn_A.sbatch` again resumes too.
Single core per run (the solver is serial FreeFem++; ~11-25 s/step locally => ~5-10 h; HPC cores should be faster). Peak RAM ~3 GB (8G requested).
Output: out_A/A_dyn_s{N}.vtu and out_F/F_dyn_s{N}.vtu every 5 steps = 0.05 ns (~300 frames each, all fields), .pvd, _log.csv, _meta.json, _ck.txt.
Restart-safe: rerunning a sbatch resumes from the last checkpoint (every 50 steps = 0.5 ns).
Post: run_post.sbatch makes 14 animations x {GIF, PNG zip, contact sheet} per variant plus raw-data zips in deliver_A, deliver_F.
Edit the module lines (freefem/4.15 sqlite/3.45.0) / python module name to match the cluster. Solver loads plugin iovtk (FF_LOADPATH is set if found).
Solver fixed vs older bundle: string literals split to compile on 4.13/4.15.

## Post-processing on the HPC (after the runs finish)
1. One-time: `bash setup_py.sh` on a login node (creates ~/dynpost_venv with numpy, scipy, matplotlib, meshio, pillow).
2. `sbatch run_post.sbatch` (or it is queued automatically by submit_all.sh). Needs out_A/ and/or out_F/ (the *_s*.vtu frames + .pvd) in the bundle folder.
3. Results: anim_{A,F}/ (14 animations x GIF and png/), deliver_{A,F}/ (GIFs, contact sheets, PNG zips, data zip, SMOOTHING.md).
Single variant / test: cd post; BASE=A_dyn NP=4 python3 animate.py ../out_A ../anim_A
Cost: ~300 frames x 14 plots; a few minutes per variant on 4 cores. Memory ~2-4 GB per worker.

## Compression + shear test (default: 10 GPa compression + 1 GPa shear, variant A, resumes from out_A/A_dyn_relaxed_ck.txt)
`bash submit_comp.sh` (chain of 6 restart-safe jobs). Env overrides: AMP, TEND, TSCALE, TAG. Output (auto tags): out_cs10_1/cs10_1_*.vtu for 10 GPa compression + 1 GPa shear, out_comp10/comp10_*.vtu for SHEAR=0 (pure compression); the original shear run stays in out_A/A_dyn_* (every 0.05 ns).
Load: traction (-AMP, +SHEAR) GPa on the right edge (label 3) with the left edge clamped (new solver option -dir combo -shear; SHEAR=0 gives pure compression); envelope = 2 ns smooth ramp, hold to 5 ns, 2 ns unload (TSCALE stretches/compresses it).
Watch min_phi / max_rhoG / mean_rhoS in the *_log.csv in the case folder for nucleation (see notes).

## Load ladder (compression P + shear 0.1 P; all levels in parallel)
`bash submit_ladder.sh [N]` -> 7 array tasks: P = 100, 200, 500, 1000, 2000, 5000, 10000 MPa (shear 10 ... 1000 MPa), each from out_A/A_dyn_relaxed_ck.txt.
Outputs out_lad_P<MPa>MPa_S<MPa>MPa/. Chain of N arrays gives restart after 12 h limits. A level that blows up (NaN) writes out_<tag>/STALLED with the end of the solver log and stops; other levels are unaffected.
Compare levels with min_phi, max_rhoG_um2, mean_rhoS_um2 in each log (nucleation shows as phi dropping where rhoS is high; threshold rho ~ 3.8e3 um^-2).
Post-process all levels: `sbatch run_post_ladder.sbatch` (array, one task per level; outputs deliver_<tag>/ with GIFs, contact sheets, PNG zips, data zip; needs setup_py.sh once).
Post-process one level: cd post; BASE=lad_P1000MPa_S100MPa NP=4 python3 animate.py ../out_lad_P1000MPa_S100MPa ../anim_lad_P1000MPa_S100MPa

## Nucleation probe
`cd post; python3 nucleus_probe.py ../out_lad_P10000MPa_S1000MPa lad_P10000MPa_S1000MPa ../nuc_P10000 [dphi_thr=0.7] [dtheta_deg=2] [frame_step=1]`
-> nuc_P10000_nucleus.csv (per frame: area of new phi-drop regions where phi was ~1 at the reference frame, area of lattice-orientation pockets
(misrotation vs grain median developed since the first frame), largest region size in nm, position, mean T and rhoS) and nuc_P10000_nucleus_last.png.
Baseline: the 15 ns shear run gives zero new phi-drop and zero pocket area (negative control).

## Confined piston ladder (cold-spray analogue), 30 ns, variants A and F  -- NEW solver option -dir piston
Set-up: wall u_x=0 (left), rollers u_y=0 (top and bottom, exactly 0), right edge prescribed u_x = -d(t); d(t) = integral of u_p*e(t), e = smooth quintic ramp up (-ptr, default 0.5 ns), hold (-pth, default 2 ns), ramp down;
piston travel = u_p*(ptr+pth) (u_p=250 m/s -> 0.625 um = 6.25 %); afterwards the piston stays put and the block evolves to 30 ns. Arguments: -pistonup (um/ns = 1000 m/s), -ptr, -pth, -tstop (K).
`bash submit_piston.sh [N]`: 12 array tasks = variants A (tasks 0-5, out_pist_up<m/s>ms) and F (tasks 6-11, out_pistF_up<m/s>ms, no relaxation stage, starts fresh), u_p = 25, 50, 100, 150, 250, 400 m/s (impact 50 ... 800 m/s; Cu shock-fit pressure 0.9 ... 16 GPa), TEND=30 ns, dt=0.01 ns (3000 loaded steps),
frames every 0.1 ns (VTKEVERY=10 steps; ~300 frames, ~13 MB each), 24 h per job, N chained arrays (default 6) resume from checkpoints. Outputs out_pist_up<m/s>ms/.
Safeguard: -tstop 1400 K (override with TSTOP=...): a level that exceeds it (or gives NaN) writes a final frame (..._tstop<step>.vtu), a checkpoint and <tag>_TSTOP.txt and ends cleanly (counted as finished).
Notes: W_ext / audit_residual columns are not meaningful in piston mode (no traction work is accumulated). The relaxed state is reused (left wall u_y freed, rollers set u_y=0 at the first piston step).
Post: `sbatch run_post_piston.sbatch` (needs setup_py.sh once); nucleus_probe.py works on these folders too (BASE=pist_up250ms).

Variant F only: `bash submit_piston_F.sh [N]` (6 array tasks, out_pistF_up<m/s>ms/) and `sbatch run_post_piston_F.sbatch` for its post-processing.

## Slip directions over lattice orientation (post-processing)
`animate.py` now also writes `stress_frob/` (local stress, Frobenius norm of sig_11, sig_22, sig_12, sig_21, sig_33 in MPa) and `lattice_slip/` (PNG frames + `lattice_slip.gif`), produced automatically by `run_post*.sbatch` and packaged by `package.sh`.
Background = `theta_lattice_deg` (period 180 deg). Overlaid short lines = the 3 slip directions at lattice angle + {0, +54.74, -54.74} deg (system 1 red, 2 blue, 3 green); line width grows with that system's |slip| (p1/p2/p3, scaled to the 99.5th percentile over sampled frames). Lines are skipped in grain-boundary cores (phi < 0.9).

## GND/recovery variant of the piston ladder (2026-10-29)

`cosserat2d_gnd.edp` is `cosserat2d.edp` with these changes and nothing else (all options default to the original behaviour):

| Change | Option | Solver default | `run_piston_gnd.sh` default |
|---|---|---|---|
| GND density without the φ² weight: ρ_G = \|∇(θ − θ_R)\|/b also in the boundary cores (original: φ²\|∇(θ − θ_R)\|/b) | `-gndphi2 0` | 1 | 0 (env `GNDPHI2`) |
| SSD recovery by boundary sweeping. 0: −C_D\|φ̇\|ρ_S (original). 1: −C_D max(φ̇,0)ρ_S. 2: −C_D max(φ̇,0) tanh(C_A²\|∇θ_t\|²) ρ_S, Ask et al. 2019 Eqs. (80), (133): recovery only behind a sweeping boundary, localised to orientation gradients | `-recov 0/1/2`, `-ca` (µm), `-cd` | 0, √10, 10 | 2, √10, 100 = Ask et al. Table 3 (env `RECOV`, `CA`, `CD`) |
| Per-step worst elements (largest slip increment, local residual, dissipation): position, grain, φ, ρ_S, ρ_G, τ*, drag | `-diag 1` → `<tag>_diag.csv` | 0 | 1 (env `DIAG`) |
| Log column `max_abs_dtheta_deg` fixed: it was max\|u\| × 180/π because `th[]` is the whole [u1,u2,θ] vector; it is now the lattice rotation θ | always | | |

The output tag encodes the switches, so nothing overwrites the original `out_pist_up*`:
`pistGR` (both, default), `pistG` (`RECOV=0`), `pistR` (`GNDPHI2=1`). Variant A only.

```
bash submit_piston_gnd.sh                    # both changes: out_pistGR_up<m/s>ms/
RECOV=0 bash submit_piston_gnd.sh            # GND only:     out_pistG_up*
GNDPHI2=1 bash submit_piston_gnd.sh          # recovery only: out_pistR_up*
sbatch run_post_piston_gnd.sbatch            # post-processing; export the same GNDPHI2/RECOV as the run
```
Single level: `bash run_piston_gnd.sh 400 [tend] [ramp] [hold] [vtkevery] A`. All need `out_A/A_dyn_relaxed_ck.txt` (same relaxed start as the original ladder).

Why: every original A piston level failed at the grain 5/grain 1 boundary where it meets the left (Dirichlet) wall, when the reflected compression (~22 GPa) builds there (2.6–3.4 ns at 400 m/s).
- **The core is weak.** φ ≈ 0.03 there, so the φ² weight removes the GND hardening.
- **ρ_S is drained.** The |φ̇| recovery term sweeps ρ_S away whenever φ changes in either direction.
- **The junction localises.** It is slip localisation with thermal feedback, and it occurs at Δt = 0.005 ns as well. The local Newton stalls are a consequence, not the cause.

## Mobile-density drag variant of the piston ladder (2026-10-29)
`cosserat2d_mob.edp` = `cosserat2d_gnd.edp` + opt-in drag options (solver defaults reproduce the original):
- `-dragrho 0|1|2|3`: mobile density in the drag eta_d = B0 (T/T0)/(b^2 rho_m): 0 rhoS+rhoG (original), 1 rhoS only, 2 Ma-Roters algebraic (comparison only), 3 evolved mobile density (Kubin-Estrin type: d rho_M/d gamma = C2[rho_Msat(1-sqrt(rho_F/rho_Fcut)) - rho_M], capped at rhoS+rhoG); `-rhomsat 100 -c2m 10 -rhofcut 1e4 -rhommin 1` (um^-2).
- `-vc <km/s>`: relativistic drag B sqrt(1+(b tau/(v_c B))^2), i.e. |gdot| <= rho_m b v_c (0 = off).
- `-relaxfix 1` (default): applies the end-of-relaxation reset (rhoS=rho0, thR=th, rhoG=0, T=T0, V=0, E0) when restarting from the relaxed checkpoint, which is written before that reset. `-relaxfix 0` = old start.
- rho_M is checkpointed in `<ck>.rhoM` (main checkpoint format unchanged). New VTU cell fields rho_mob, eta_drag; new diag columns rho_mob, eta_eff, tau_abs_max, v_kms.
Run: `bash submit_piston_mob.sh` (defaults DRAGRHO=3 VC=2 GNDPHI2=1 RECOV=2 CD=100 -> out_pistM3v2RcV_up<m/s>ms/; since 2026-10-30 the default recovery is Ask 2019 Eq. 80, RECOV=0 gives the earlier pistM3v2 runs); post: `sbatch run_post_piston_mob.sbatch` with the same env. animate.py adds rho_mob/eta_drag animations automatically when the fields are present.

### Robustness fixes in cosserat2d_mob.edp (2026-10-30)
Diagnosis of the ~6-7 ns stops (pistM3v2_dt0.005_nr10_up25ms, step 1365): not physical heating (max dT 0.9 K before) but one non-converged global step (40 its, residual 46) committed anyway. One boundary element at the top edge took a 2.9 rad spin increment; the q^4-series exponential in trialFp() is not a rotation at that size, so det Fp != 1, J = 0.42, stresses ~110 GPa, T = 42000 K.
- `-exactexp 1` (default): exact exp of the traceless 2x2 increment, exp(A) = c I + s A with c = cosh q / cos w, s = sinh q/q / sin w/w (Cayley-Hamilton); series only for |q| < 1e-3. det Fp = 1 for any increment. Identical to the old update for small increments.
- `-spinmax <rad>` (default 0 = off): bound |dw*| per step.
- `-onfail retry|stop|commit` (default retry): a step that hits maxit (or gives NaN) is restored to its predictor, the tangent is refreshed there and the step re-solved with up to `-maxitretry` (80) iterations; if that fails too, the run stops cleanly: checkpoint at the end of the previous step, `<out>_stall<n>.vtu`, `<out>_NOCONV.txt`. `commit` = old behaviour. The log `iters` column of a retried step is the retry's count.
- `-ckroll N` (scripts: 20): rolling checkpoints `<out>_ckroll0.txt` / `_ckroll1.txt` every N steps, so a crash can be replayed (`-restart <out>_ckroll1.txt`).
Scripts: env SPINMAX (tag suffix _sm<value> when set), ONFAIL, MAXITRETRY, CKROLL, EXACTEXP. run_piston_mob.sh refuses to resume a level with a `_NOCONV.txt` until it is deleted.

### Thermodynamic consistency (2026-10-30, see clausius_duhem.pdf)
Coleman-Noll gives D = sum tau_a gdot_a + tau* w*dot + eta_phi phidot^2 - Z rhodot_S (+ conduction), Z = chi mu b^2 phi / 2. The continuum kinetics satisfy D >= 0 (storage uses <= 8% of the plastic work below ~2400 K). The negative-D counts came from the time discretization: the heat source used the midpoint stress (tau_n + tau_{n+1})/2, which is negative when the stress drops within a step, and from committed non-converged steps.
- `-dissip 1` (scripts: DISSIP=1): heat source tau_{n+1} dgamma (sign-definite for the backward-Euler local update); 0 = old midpoint.
- `-cdstore 1` (scripts: CDSTORE=1): KMT storage capped so that Z drho_store <= sum|tau dgamma| (inactive in normal flow).
- `-softmin 0.01` (always): lower bound of the thermal softening factor, keeps tau_c > 0.
- New `<tag>_thermo.csv`: per-step heat, implicit / midpoint slip dissipation, cumulative midpoint-minus-implicit, number of elements with D < 0, min D, stored-energy rate.
Tag: a `c` is appended when DISSIP or CDSTORE is on (default), e.g. out_pistM3v2RcV_dt0.005_nr10_up400ms.

### Elastic law: Vinet / Murnaghan equation of state (2026-10-30)
SVK on Ebar loses stiffness under compression (dP/dlambda = C(3 lambda^2 - 1)/2 = 0 at lambda = 1/sqrt 3; max uniaxial compressive stress ~0.19 C = 32-42 GPa) and an element pushed past that collapses to J -> 0.
`-eos vinet|murnaghan|svk` (solver default svk; scripts default EOS=vinet): psi_e = W(J) + 1/2 Et:C:Et + 2 muc a^2 with the isochoric Green strain Et = (J^{-2/3} Cbar - I)/2 and the full 3D cubic C in plane strain (C_ij33 = B_ij from Clxz, Clyz; C_3333 = Clxx).
Vinet: P = 3 K0 (1-x)/x^2 exp(eta(1-x)), x = J^{1/3}, eta = 1.5(K'-1); Murnaghan: P = K0/K' (J^-K' - 1). K0 = (C11+2C12)/3 = 137.1 GPa (from the elastic constants, so the small-strain response is unchanged), K' = 5.0 (`-kp`, script KP).
`-voltan 1` adds (K(J) - K0)+ (div v)(div w) to the tangent at every refresh (the stiffening otherwise stalls the modified Newton iteration).
Checks: S = dW/dE against finite differences (python); the solver's -check finite-difference tests pass with -eos vinet (slip driving forces to 1e-5, internal force to 1.5e-4); uniaxial-strain P(lambda) is monotone (no collapse).
Tag: V (vinet) or Mu (murnaghan) after the R/c letters, e.g. out_pistM3v2RcV_dt0.005_nr10_up400ms.

### Stage breadcrumbs and NaN guards (2026-10-31)
`-trace 1` (default): `<tag>_stage.txt` always holds the current step, global iteration and stage (step start / predictor / local plastic update / residual / linear solve / tangent: plastic condensation / tangent: assemble+factorize / commit: Newmark+dissipation / phase-field solve / SSD+mobile density / heat solve / output). After a hard crash (exec error, UMFPACK failure, kill) it shows exactly where the solver was.
NaN/Inf guards after each stage write `<tag>_NAN.txt` (stage, array, first bad element with grain, centroid, phi, J). A NaN inside the global iteration stops at once (no retry, so no factorization of a NaN matrix) with a checkpoint at the end of the previous step; a NaN during the commit stops without overwriting any checkpoint. run_piston_mob.sh prints the last stage when the solver exits with an error and refuses to resume a level with a `_NAN.txt` until it is deleted.

### Global Newton iteration (2026-10-31)
`-newton full` (scripts default NEWTON=full; solver default `modified` = old behaviour): the tangent (Jacobian) is rebuilt at the current iterate in every global iteration (KWC curvature tangent, plastic condensation at the current slip state, volumetric EOS correction), the correction is scaled by a backtracking line search on ||R|| (Armijo, halving, at most `-lsmax` 6), and a step converges when du < `-tol` (1e-7) and ||R|| <= `-tolres` (1e-6) * max(||R_0||, `-resfloor` 1e-8). If the increment becomes negligible (du < 1e-9) while the residual sits on the floor of the inexact local solve, the step is accepted and counted ("converged on du only" in the .out). `-nlocal` (30) sets the per-element slip/spin Newton iterations. With -newton full the scheduled refresh (NREF) is not used.
Local test (400 m/s, t = 2.5 ns): residual 121 -> 4.8e-6 in 11 iterations (line search used in ~1 of 5 iterations); the converged state agrees with the modified-Newton one to all printed digits; peak memory 3.95 GB (no growth over 22 refactorizations); ~6x the wall time per step of the modified iteration (~400 s vs 70 s locally). For the full 30 ns ladder use more chained arrays: run_all_A.sh now defaults to 12.
The loop is one function (newtonSolve) for the normal solve and the retry; the step logs `|R|/|R0|` in the NOCONV message and in verbose output.

### Fully consistent tangent (2026-10-31)
`-tangent consistent` (scripts default TANGENT=consistent; solver default `approx` = small-strain stiffness + small-strain plastic condensation). For each element the stress s = (P11, P12, P21, P22, dpsi/dtheta) and the driving forces tau = (tau_1..3, tau*) are closed-form functions of the element-mean kinematics k = (F11-1, F12, F21, F22-1, theta) and the local increments q = (dgamma_1..3, dw*). Their partial derivatives ds/dk, ds/dq, dtau/dk, dtau/dq are taken by central differences of that closed-form map (no local solve inside; error ~h^2, h_k = 1e-6, h_q = 1e-8); the local equations r_a = dgamma_a - dt sgn(tau_a) g_a(|tau_a|), r_w = dw* - kspin tau* are differentiated analytically (g' from the flow law incl. the relativistic drag), dq/dk = -(dr/dq)^-1 dr/dk is solved per element (4x4), and D_tot = ds/dk + ds/dq dq/dk minus the quadratic-energy element stiffness (already in the pointwise operator) is assembled as Kk' diag(D/A) Kk with the element-mean gradient operators. The KWC curvature tangent stays analytic. Covers the finite-strain elasticity, the EOS, the exact exponential, the Cosserat skew, slip with drag and velocity cap, and the coupled GB spin.
The same exact 4x4 local Jacobian is used for a coupled slip + spin Newton (`-polish`, 20 iterations max, per-element backtracking) after the robust local loop, so the local state is converged to round-off.
Checks (400 m/s, t = 2.5 ns): directional finite-difference check of the global residual (`-tancheck 1`): relative error 4.2e-5 with the consistent tangent vs 7.1e-2 with the approximate one. Terminal convergence quadratic: |R| 6.9e-4 -> 1.3e-5 -> 9.9e-9 -> 1.8e-11 (|R0| = 121). The first iterations are governed by the m = 50 rate-law nonlinearity and use the line search.
Cost: ~10x the modified Newton per step; dominated by the assembly and UMFPACK factorization of the global matrix every iteration (the finite-difference element derivatives are ~1 s per tangent).

### Speed-ups (2026-10-31)
- `-linsolver gmres` (scripts default LINSOLVER=gmres): restarted GMRES (Krylov 100, relative tolerance `-lintol` 1e-10) with an ILU(0) preconditioner (IncompleteCholesky plugin). The ILU is rebuilt once per step (`-precevery 0`) and reused between tangent rebuilds. Every solve's residual is checked; a solve that misses 1e3*lintol is redone with the direct solver (counted as "GMRES fallbacks"). run_piston_mob.sh falls back to umfpack automatically if the plugin is missing.
- The constant part of the matrix is assembled once; with the consistent tangent the whole matrix is assembled in one varf per rebuild (no sparse matrix additions; the element tangent uses the barycentre quadrature, which gives the element-mean gradients exactly).
- `-tanevery k` (default 1) rebuilds the tangent every k iterations; `-slippred 1` starts the local solve from the previous step's increments. Neither helped on the test step (k = 2, 3: 14-18 iterations instead of 11, ~same wall time; slippred: same iteration count), so both are off by default.
Measured, one 400 m/s step at t = 2.5 ns, full Newton + consistent tangent + Vinet, 11 iterations: 521 s -> 120 s (linear solve 233 -> 36 s, tangent/assembly 171 -> 56 s); same residual history and converged state to all printed digits. Old modified-Newton path unchanged (reproduces the earlier runs).
Checked end to end with run_piston_mob.sh on 25 m/s from the relaxed checkpoint (relaxfix, Vinet, full Newton, consistent tangent, GMRES; 3 iterations per step in the elastic start), including a resume with log trimming.

### NaN diagnosis (2026-10-31) and guards
All three NaN stops (50, 250, 400 m/s; NREF=1 runs) are the same event: J = det Ubar becomes NaN during the local plastic update in global iteration 4-6 of the step at t = 2.79-2.85 ns, in grain 8 next to the top boundary (y = 9.66-9.96 um). The slip/spin increments checked just before were finite, so the global iteration had already diverged and the huge iterate overflowed cosh/sinh in the exact exponential. Guards added: the exponential argument is capped at 50 (finite values; such an iterate is then rejected by its residual), and an iteration whose residual exceeds 1e4 |R0| is abandoned at once (-> retry with a fresh tangent / clean stop) instead of being iterated into NaN.
probe_mpi.edp / probe_mumps.edp: check whether the cluster FreeFem has PETSc (route for the MPI version) or parallel MUMPS.

## MPI/PETSc version (2026-10-31): cosserat2d_mob_mpi.edp
Same physics, options, defaults and output formats as cosserat2d_mob.edp. It is generated from it, so the serial solver stays the reference:
`python3 mpi/make_mpi.py cosserat2d_mob.edp cosserat2d_mob_mpi.edp` (rerun this after any change to the serial file; it stops with an error if a patched line has changed).
* Domain decomposition (FreeFem macro_ddm / buildDmesh, overlap 1): each rank does the element work (local slip/spin Newton, consistent tangent, assembly) on its own subdomain.
  The (u, theta) Newton system, the phase field and the heat equation are PETSc solves. Norms, energies, maxima and every convergence/termination test are global reductions, so all ranks make the same decisions.
* Rank 0 writes the logs, diag, thermo, VTU and checkpoints on the global mesh in the serial formats: post-processing is unchanged, and checkpoints are interchangeable both ways between serial and MPI.
* Verified on step 351 of the 400 m/s run (from a 2.5 ns checkpoint, full Newton + consistent tangent):

  | | serial | 2 ranks | 4 ranks |
  |---|---|---|---|
  | Newton history | 11 its, 1.84e-11/142.765 | identical | identical |
  | checkpoint difference vs serial | – | 2e-9 | 2e-9 |
  | VTU fields | – | equal to float32 precision | – |
  | time per step (CPU per rank) | ~220 s | 88 s | ~50 s |
* Solver options. `-petsc` (default `-ksp_type gmres -ksp_gmres_restart 200 -ksp_rtol 1e-10 -pc_type asm -sub_pc_type ilu`, i.e. additive Schwarz with ILU(0) per subdomain) is used for the (u, theta) Newton system. Subdomain LU was 14x slower in the solve; MUMPS direct was also slower.
  `-petscfall` (default MUMPS direct) redoes any solve that misses -lintol. `-petscp1` sets the solver for the phase-field and heat systems. `-linsolver` is ignored.
* Scripts:
  * run_piston_mob_mpi.sh / run_piston_mob_mpi.sbatch take the same environment options as the serial scripts, plus:
    * `NP` ranks (sbatch default 8, 24 GB);
    * `MPIRUN` (default: the mpiexec recorded in `ff-mpirun`, i.e. the MPI FreeFem++-mpi was built with; else `mpiexec`/`mpirun`, else `srun`; the solver stops at once with LAUNCHER ERROR if it does not get NP ranks);
    * `PETSC`, `PETSCFALL`, `PETSCP1`;
    * `SHARETAG=1`, which resumes the serial run of the same tag instead of using the separate `_mpi` tag.
  * run_all_A_mpi.sh chains the six speeds (default DT=0.01, NREF=20, 10 links) plus run_post_piston_mob_mpi.sbatch:

        RUNTAG=gm bash run_all_A_mpi.sh            # -> out_pistM3v2RcV_gm_mpi_up<speed>ms/
        RUNTAG=gm DT=0.005 NREF=10 NP=16 bash run_all_A_mpi.sh
* Needs `FreeFem++-mpi` with the PETSc plugin (the freefem/4.15 module: check with `mpiexec -n 2 FreeFem++-mpi -nw probe_mpi.edp`) and the MPI launcher it was built with on PATH.

## GB-spin branch fix (2026-10-31): -spinbranch / SPINBRANCH
Symptom: a single grain-boundary element at the free top surface had its GB spin jump in one converged step. In the
50 m/s MPI run at 3.96 ns: element 1068, grain 8, phi = 0.6, w* went from 0.001 to 71.8 rad, and its neighbour to -3.08 ≈ -pi.
The heat source there exploded and the run hit TSTOP at T = 6e5 K. This is the same location and speed independence as the earlier serial NaN stops.
It reproduces identically on 2 and 8 ranks, so it is not an MPI effect.
Cause: with the exact exponential, the spin increment is an exact rotation of the plastic frame. The energy is pi-periodic in w*
(Ubar -> -Ubar leaves Ubar^T Ubar and the Cosserat term a^2 unchanged), so the local equation dw* = k tau*(dw*) has spurious roots near n*pi
wherever k = dt (1-phi)^2/(2 eta*_0) is large (boundary cores, k ≈ 0.8). The local Newton at the predictor, and then the polish, landed on one.
Fix: `-spinbranch b` keeps |dw*| <= b per step in both the spin update and the coupled polish (solver default 0 = old behaviour;
script default SPINBRANCH = 0.785398 = pi/4). Physical increments are ~1e-3 rad per step (the largest accumulated |w*| in that run is 0.023 rad), so the bound never binds a converged state.
It only keeps the intermediate iterates on the principal root. Replay of steps 481-500:
* steps 481-495 are identical with and without the bound;
* step 496 converges in 13 iterations with max |dw*| = 3e-4 rad and max dT = 2.6 K;
* steps 497-500 are normal.
Not in the run tag (it does not change any converged physical state), so stopped runs resume under the same tag:

    bash recover_ck.sh out_<tag> [...]     # rolls <tag>_ck.txt back to the newest rolling checkpoint (the TSTOP/NaN stop saved the bad state), moves the stop markers to *.old
    # then resubmit as before; SPINBRANCH=0 restores the old behaviour

## Prescribed-dof kinematics and the Newton predictor (2026-10-31): -bckin / BCKIN, -predictor / PREDICTOR
Symptom: the 50 m/s MPI run stopped (NOCONV) at 6.23 ns. The global Newton diverged from its very first iterate:
the predictor already put J = 0.49 next to the piston face, and the first correction inverted elements (J = -1088).
Cause: on the prescribed dofs (piston face u_x, rollers u_y) the Newmark update a_{n+1} = c0 (U-Un) - c1 v_n - a_n was still applied.
With U - Un = 0 it makes v flip sign every step and a grow linearly at 4|v|/dt per step (measured 7.49/step against 4*0.0187/0.01 = 7.48).
The kick comes from the rollers being snapped to 0 after the relaxation. Through the consistent mass these dofs push their free neighbours with
a growing alternating force: at 6.2 ns |a| = 4.7e3 um/ns^2 on the prescribed dofs and 1.2e3 on their free neighbours.
The average-acceleration predictor turns that into a 0.03 um kink (dt^2/4 a) before Newton starts.
This is likely also the trigger of the 3.96 ns GB-spin blow-up (an element on the top roller boundary, where the predictor's local residual was already 365).
It probably also explains the earlier "crash at ~7 ns" runs.
Fixes (solver defaults = old behaviour; script defaults BCKIN=1, PREDICTOR=vel):
* `-bckin 1`: after every commit, and once on restart (cleans existing checkpoints), the prescribed dofs carry their exact kinematics:
  the piston face gets -v_p(t), -a_p(t) (derivatives of the quintic ramp-and-hold); rollers and the wall get 0.
* `-predictor vel`: Newton starts from U + dt V instead of U + dt V + dt^2/4 A. This is only the starting guess: with full Newton and the consistent tangent the converged step does not depend on it.
Replay from the step-722 checkpoint (2 ranks):
* both off: diverges at step 723;
* `-bckin 1` alone: still diverges, because the free-dof acceleration left over from the spurious mode still builds the kink;
* `-predictor vel` alone: converges, but the spurious boundary acceleration keeps growing (4685 at step 726);
* both: steps 723-730 converge in 6-10 iterations each; the prescribed |a| is 0 and the free-dof residue is no longer driven.
Not in the run tag, so stopped runs resume under the same tag (resume_all_A.sh clears the NOCONV marker; the NOCONV checkpoint is clean).
Serial copies already patched with the spin fix can be updated in place with `bash apply_bckin_serial.sh` (atomic, safe while jobs run).

## Angled plate impact (2026-11-01): -dir shearpiston, ANGLE
The same plate-speed ladder V = 25 50 100 150 250 400 m/s, with the plate velocity at an angle ANGLE (degrees) from -x:
u_n = V cos(ANGLE) (compression, -x) and u_t = V sin(ANGLE) (shear, +y). ANGLE = 0 is pure compression, 90 is pure shear; the script default is 45.
Boundary conditions (the same for every angle, so the angle is the only variable):
* left edge clamped (u = 0, no-slip substrate);
* right edge the plate: u_x = -d_n(t), u_y = +d_t(t), same quintic ramp / hold / ramp-down as the piston (RAMP, HOLD), exact plate kinematics via -bckin;
* top and bottom edges traction-free.
ANGLE = 0 here is therefore NOT the old confined piston (that had rollers on the top and bottom edges and a frictionless wall), and the two are not comparable one to one.
ANGLE=confined runs the old confined piston under the old tag (no _a), so existing runs resume.
Solver options: -dir shearpiston -pistonup <u_n um/ns> -pistonut <u_t um/ns>. The run scripts compute both from V and ANGLE. Tag: ..._a<ANGLE>_up<V>ms.
Notes:
* 400 m/s of shear needs a plastic rate of ~1e8-1e9 /s in a band. The relativistic-drag ceiling rho_m b v_c is ~5e7 /s per system at the default
  RHOMSAT = 100 um^-2 (1e14 m^-2), so the response at high shear speeds is largely elastic unless RHOMSAT is raised (1e15-1e16 m^-2 = 1000-10000 um^-2 is typical of high-rate Cu).
* The tangential travel is V sin(ANGLE) (RAMP + HOLD) (1 um at 90 deg / 400 m/s); localized, that is a shear strain of O(1).
  The Lagrangian mesh distorts there, and a shorter HOLD is the cautious first run.
Examples (MPI):

    RUNTAG=bk bash run_all_A_mpi.sh                 # 45 deg  -> out_pistM3v2RcV_bk_mpi_a45_up<V>ms
    RUNTAG=bk ANGLE=90 bash run_all_A_mpi.sh        # pure shear
    RUNTAG=bk ANGLE=0 bash run_all_A_mpi.sh         # pure compression with the plate BCs
    MPI=1 RUNTAG=bk ANGLE=45 bash resume_all_A.sh   # resume; ANGLE=confined for the earlier confined-piston bk series

### tau / temperature / mismatch statistics in the post step (added 2026-11-08)
`run_post_piston_mob.sbatch` and `run_post_piston_mob_mpi.sbatch` now also run `post/tau_stats.py` on `out_<tag>` (after grain_post) and write `tau_<tag>/` (per-snapshot CSV, `tau_summary.csv`, `tau_vs_time.png`); `package.sh` copies it to `deliver_<tag>/tau/` (optional 6th argument). It reports tau_1..3 and tau_c, T, pressure, the heating-jump time/rate, the Cosserat mismatch beta = phi_e - theta and the solver log's `max_abs_dtheta_deg`. Env: `TAUWIN0/TAUWIN1` (loading window, default 0.5-30 ns), `TAUMARK` (event marker, default 6 ns), `TAUEVERY` (every N-th snapshot). A failure does not affect the animations. Only syntax-checked (sbatch) and tested on synthetic data (package.sh, tau_stats.py); not run on Longleaf.
