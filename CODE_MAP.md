# Code map

Where things are in the code. How to run it: `README.md`. What each output column and field means: `DATA_FIELDS.md`.

Version 8.0 is the production freeze: the files listed in `SOLVER_SHA256` are not edited any more (`checks/verify_solver.sh`).

## Solver `cosserat_kwc.edp`

The single source of the solver. Sections in file order; search for the function name (line numbers drift with edits).
The MPI file `cosserat_kwc_mpi.edp` is generated from it by `mpi/make_mpi.py` (below); never edit the generated file.

| Section | Key functions and variables | Contents |
|---|---|---|
| header | | summary, units, plugins: `iovtk`, `IncompleteCholesky`, `kwc2d` (this repository's plugin, `plugin/`) |
| command-line arguments | `arg()` | run options (`-pistonup/-pistonut/-ptr/-pth/-tend/-dt/-out/-restart/-walllimit/…`) |
| material parameters | `gmaxS`, `gsumS`, `gdotS`, `isRoot`, `syncPen` | elastic constants, Vinet `-kp/-k0`, slip systems, KWC `-f0/-kwca/-kwcs/-kwcd`, φ kinetics, recovery, drag, mobile density; serial stand-ins of the MPI collectives |
| mesh and grain data | `Wh`, `Vh1`, `Ph` | `mesh.msh`, `grain_angles_rad.txt`, P2 displacement + P2 θ space `Wh`, P1 φ space, P0 history space, plate and clamp dofs; per-element constants handed to the plugin (`k2dSetup`, `k2dThreads`) |
| history, initial orientation | | P0 arrays (slip, GB spin, densities, F_p); smooth initial θ field |
| operators | `buildA`, `linSolve` | mass matrices, KWC operators, GMRES + ILU with a direct fallback |
| Gauss-point kinematics | `jBound` | `D1..D5` (U → ∇u, θ at the 3 edge midpoints), det F Bernstein bound (`-jmin`) |
| local constitutive update | `effDrag`, `flow`, `dragBound` | relativistic drag, flow rule, small solvers |
| finite-deformation stresses | `trialFp`, `eosStress`, `finiteStress`, `taus`, `backtrack`, `kin` | trial F_p, Vinet + isochoric energy, stresses, resolved stresses, local backtracking, centroid kinematics |
| consistent tangent and local Newton | `solve4`, `elemJac`, `buildJqq`, `localRes`, `polishLocal`, `localUpdate`, `strength` | `elemJac`: derivatives of the element map (stresses and resolved stresses w.r.t. F, θ and the slip/spin increments) by forward-mode AD in the plugin (`k2dJac`); local slip + spin update with the exponential map and the \|q\| cap |
| state vectors, plate | `fixBCkin`, `kwcStress` | plate travel/velocity/acceleration (quintic ramp, hold), φ² weights, KWC couple stress |
| internal force | `gaussInt`, `residual` | f_int = dE/dU at the Gauss points, energy `gaussE` |
| exact Gauss Hessian | `packPd`, `gaussTangent`, `consistentTangent` | one-form `vJg` assembly of the tangent including the implicit slip/spin response |
| guards and energies | `mark`, `nanE`, `nanV`, `refreshTangent`, `elasticEnergy`, `storedEnergy` | stage breadcrumbs (`<out>_stage.txt`, `-trace`), NaN guards, kinetic energies, heat |
| output | `writeVTK` | P0 → P1 projection, VTU + PVD |
| GND | `curlFpG` | Cermelli–Gurtin, planar: g = c/J_p, ρ_G = \|g\|/b |
| time-loop setup | `mobUpdate`, `linSolve` | step counts (`nrelax + tend/dt`), mobile-density update |
| global Newton | `newtonSolve`, `plateLog`, `reactionNow` | full Newton, line search on \|R\|, stagnation exit; plate work ledger |
| checkpoints and restart | `saveCkF`, `saveCkA` | atomic generations with manifest and `.gnd/.plate/.rhoM` side files; restart reader (a relaxed checkpoint without side files is accepted) |
| run metadata | | `<out>_meta.json`, CSV headers |
| time loop | | per step: φ weights and drag → predictor (U + dt V, plate; exact plate acceleration in the inertial load of the prescribed dofs) → Newton (retry with a fresh tangent, `_NOCONV.txt`) → Newmark commit and plate work → dissipation → φ solve (primal-dual active set) and its dissipation (element integral of η_φ(T)(φ − φ_n)²/Δt, `vDph`) → SSD storage/recovery, GND, mobile density → heat. At `n == nrelax`: relaxed checkpoint `<out>_relaxed_ck.txt` and `_relaxed.vtu`, then the end-of-relaxation reset (V, A, T, ρ_S, θ_R, G_R). Output: log/thermo/diag rows, VTU frames, checkpoints, rolling checkpoints, TSTOP and wall-limit stops |

## Plugin `plugin/`

- `kwc2d.cpp`: the element map of the solver (same formulas as `trialFp`, `eosStress`, `finiteStress`, `taus`) and its exact derivatives
  by forward-mode automatic differentiation, OpenMP over elements. FreeFem functions: `k2dSetup` (per-element constants), `k2dThreads`,
  `k2dJac` (derivatives, the solver's `DSF`/`DTQ` layout), `k2dMap` (values only; used by the `-adcheck` self-test), `k2dLam` (Cosserat
  multiplier per element; v2 solver only, without it the element map is unchanged), `k2dJacL` (as `k2dJac` with the
  multiplier as a tenth differentiated input and the mismatch a as an extra output; `-almc 4` of the v2 solver), `k2dAsmPattern/Setup/J/A/Probe/ProbeAbs`
  (direct assembly of the tangent matrix on a fixed pattern, `-asm 1` of the v2 solver).
- `dual.hpp`: the dual-number type.
- `build_plugin.sh`: builds `kwc2d.so` into the repository root (and a copy in `mpi/`, where FreeFem++-mpi looks first) and regenerates
  `cosserat_kwc_mpi.edp` when stale. `run.sh` and the submit scripts call it; it rebuilds only what changed.

A change to the element map in the solver must be made in `kwc2d.cpp` as well; `checks/run_tangent_check.sh` (-adcheck) catches a mismatch.

## Production solver `v2/make_v2.py`

Exact text patches on `cosserat_kwc.edp` giving `cosserat_kwc_v2.edp` (then `mpi/make_mpi.py` for its MPI version). With its defaults
(`-kwctan 0 -almc 0 -asm 0`) it takes the same steps as the base solver; production runs (both set-ups) pass `-kwctan 3 -almc 4 -asm 1` and `-slipexact 1 -localtol 1e-10 -betamax 0.785398`. Where the patches act:
- declarations after the KWC tangent fields: options, dual flux w, multiplier λ, constraint blocks (`DAL`, `DQL`, `cBm`, `cCm`, `eCon`);
- `finiteStress`: λ a in ψ and λ/2 in the skew stress; plugin call `k2dJacL` instead of `k2dJac` when `-almc 4`;
- basis gradients at the Gauss points and the centroid (`gbX`, `gbY`) and the element data for the fast assembly (`geoA`, `connA`);
- `consistentTangent`: sixth right-hand side dq/dλ; `gaussTangent`: coefficients of B, C, E; with `-asm 1` the plugin's `k2dAsmJ`;
- `refreshTangent`: primal–dual KWC tangent (`-kwctan 2` symmetrised; `-kwctan 3` adds the antisymmetric part `Hd`, also in the form `vKv` and in `k2dAsmA`, fourth column); `buildA`: `k2dAsmA` on the fixed pattern and the z'Ax check (serial: the matrix is rebuilt
  from its triplets once a direct solver has been attached to it, `asmSerialDirect`);
- `powerSplit` (accepted state: plastic power with centroid driving forces against their Gauss-point mean; end-of-run line);
- `flowExact` (`-slipexact`), `localAcceptRes` (`-localtol`); `betaNow` (max |β| of the current Ū arrays): in `gaussInt` at the three Gauss points, in `residual` with the element means (`-betamax`, inadmissible iterates); `acceptFinal`: the final state of a converged step (one local update, residual, constraint, local and branch checks), committed without a further local update (`stateFinal`);
- `conRes`, `linSolveLM` (saddle-point solve; replaced by the PETSc block solve in the MPI file, markers `//LM-SOLVE-BEGIN/END`);
- `newtonSolve`: merit with the constraint, U and λ in the line search, dual-flux update, convergence on max |ā|, `<out>_newton.csv`
  (max |ā| and max |β| per step); local and branch acceptance tests after convergence, in the first attempt and the retry;
- time loop: w and λ saved at the step start and restored on a retry, whole-step timing line;
- checkpoints: side file `.mult` (λ, w), written and read; metadata entry `"v2"`.

## Fast assembly of the staggered fields and the phase solve (version 7.6, in `cosserat_kwc.edp`)

- options `-fastasm` (`fastAsm`), `-fastcheck`, `-fastchecktol`; `fastDue` (which calls are compared), `phRegime` (first call with the
  Arrhenius floor / a node on a bound), `relRow` (row-wise relative comparison), `absMat`, `fastStop` (`<out>_FASTCHECK.txt`);
- phase field: plugin `k2dPhSetup` (pattern), `k2dPhAsm` (matrix and residual, 7-point rule), `k2dPhElem` (element means of φ², g(φ),
  phase dissipation); diagonal $A\mathbf 1$ and right-hand side $A\phi_k-R$ by algebra;
- GND: `CurlA`, `CurlB` (constant P0 → P1 operators), `curlFast`; heat right-hand side by constant operators;
- phase linear solve: `-phsolver cg` (`phCG`), `phPCG` (CG + Jacobi on the eliminated system, direct fallback, `nPhCGfail`), a failed
  phase step repeated once with the previous solver (`nPhFall`); MPI: PETSc CG + Jacobi, option `-petscph` (`mpi/make_mpi.py`);
- storage cap: `nStoreCapEl`, `nStoreCapSteps`, `EstoreCapCum`, `<out>_storecap.csv`.

## MPI generator `mpi/make_mpi.py`

Exact text patches on the serial file (each must match the stated number of times):
PETSc `macro_ddm` setup and the overlapping decomposition; the collective wrappers (`gmaxS`, `gsumS`, `gdotS`, `gsumEs`, `isRoot`,
`syncPen`, `syncP1`) replaced by `mpiAllReduce` versions; owned-element masks in every element sum (energies, dissipation, diagnostics);
the linear solves through PETSc (`set(A, sparams = …)`; options `-petsc` (u–θ system: GMRES + ASM/ILU), `-petscfall` (MUMPS fallback),
`-petscp1` (φ and heat), set from the environment variables `PETSC`, `PETSCFALL`, `PETSCP1` by `run.sh`); the rank-count guard (anchored
on the `-j` option line); root-only output and checkpoint writes with gathered P0/P1 arrays; `"mpi"` entry in the metadata.

## Scripts (call graph; what each does and how to use it: README.md)

```
submit_ladder.sh ─┬─ plugin/build_plugin.sh
                  ├─ run_relax.sbatch ── run_relax.sh ── run.sh (RELAXONLY=1) ── tag.sh
                  ├─ submit_speeds.sh (per angle) ─┬─ plugin/build_plugin.sh
                  │                                ├─ run.sbatch × N (array; afterok relax, then afterany) ── run.sh ── tag.sh, plugin/build_plugin.sh
                  │                                └─ post.sbatch (array) ── tag.sh, post/animate.py, grain_post.py, tau_stats.py, package.sh
                  └─ grain_summary.sbatch ── tag.sh, post/grain_summary.py
post_all.sh ─┬─ post.sbatch (array, per angle)
             └─ grain_summary.sbatch (afterany all post arrays)
resume.sh ── recover_ck.sh, then submit_speeds.sh
```

MODE (mpi | serial, default mpi) is read by every script and passed on through `sbatch --export=ALL`.

## Tools `tools/` (not used by the runs)

- `kwc_fit/`: `kwc_fit_tschopp.py`, `tschopp_grab.py`, `kwc1d.py`, `README_tschopp.md` (fit of the KWC constants to Tschopp et al. 2015).
- `debug/`: `make_dbg.py` (debug build of the MPI solver, `cosserat_kwc_v2_dbg_mpi.edp`); `diag_step.sh` (repository root) reruns
  steps after a checkpoint with per-iteration output.

## Checks `checks/`

- `make_checks.py`: builds `cosserat_kwc_check.edp` (repository root) = solver + three self-tests: `-adcheck 1` (AD element-map derivatives
  against central differences on every element, plugin values against the solver's; the central-difference element map lives only
  here), `-tancheck 1|2` (assembled tangent against central differences of the global residual, FD-step scan, KWC tangent alone;
  inserted before `newtonSolve`), `-check 1` (driving forces and internal force against energy derivatives, appended).
- `run_tangent_check.sh`: runs all three, serially, on the v2 solver with the exact KWC Hessian, from `relax_checkpoint_v2.txt` (`OLD=1`: `relax_checkpoint_j.txt`).
- `compare_runs.py`: per-step comparison of two runs of the same set-up (log columns, Newton iterations, last common VTU, wall times).
- `probe_mpi.edp`, `probe_mumps.edp`: installation probes for a FreeFem module.
- `verify_solver.sh`: `sha256sum -c SOLVER_SHA256` (the frozen solver files).
- `phase_solver/`: `kkt2.py` (Python port of the phase step, KKT solver comparison), `bench.edp` (FreeFEM timing of the linear solves),
  `ind.py` (independence analysis of the periodic matrix check); `gnd_cg_test.edp`: test of the Cermelli–Gurtin GND measure.

## Installation `cluster/`

Building FreeFEM with PETSc as a module on Longleaf and Sycamore, and testing it: `cluster/README.md`.
