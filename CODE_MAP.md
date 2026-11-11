# Code map

What each output field means: `DATA_FIELDS.md`.

Line numbers refer to `cosserat_kwc.edp` (1639 lines) as of the cleanup; they drift with edits, so search for the section header.
The MPI file is generated from this one by `mpi/make_mpi.py`; the generator's patches are listed at the end.

## Solver `cosserat_kwc.edp`

| Lines | Section | Contents |
|---|---|---|
| 1–9 | header | summary, units, plugins (`iovtk`, `IncompleteCholesky`) |
| 10–43 | command-line arguments | `arg()`, run options (`-pistonup/-pistonut/-ptr/-pth/-tend/-dt/-out/-restart/-walllimit/…`) |
| 44–98 | material parameters | elastic constants, Vinet `-kp/-k0`, slip systems, KWC `-f0/-kwca/-kwcs/-kwcd`, φ kinetics, recovery, drag, mobile density; serial stand-ins of the MPI collectives (`gmaxS`, `gsumS`, `gdotS`, `isRoot`, `syncPen`) |
| 99–170 | mesh and grain data | `mesh.msh`, `grain_angles_rad.txt`, P2 displacement + P2 θ space `Wh`, P1 φ space, P0 history space, plate and clamp dofs |
| 171–194 | history and initial orientation | P0 arrays (slip, GB spin, densities, F_p), smooth initial θ field |
| 195–267 | operators | mass matrices, KWC operators, `buildA()`, GMRES + ILU (`PreconILU`) with direct fallback |
| 268–346 | Gauss-point kinematics | `D1..D5` (U → ∇u, θ at the 3 edge midpoints), det F Bernstein bound `jBound` (`-jmin`) |
| 347–423 | local constitutive update | relativistic drag `effDrag`, flow rule `flow`, `dragBound`, small solvers |
| 424–594 | finite-deformation stresses | `trialFp`, Vinet + isochoric energy `eosStress`, `finiteStress`, resolved stresses `taus`, local `backtrack`, centroid kinematics `kin` |
| 595–770 | consistent tangent and local Newton | `solve4`, `fdMap`, `buildJqq`, `localRes`, `polishLocal`, `localUpdate` (slip + spin, exponential map with the \|q\| cap), `strength` |
| 771–842 | state vectors, plate | plate travel/velocity/acceleration (quintic ramp, hold), `fixBCkin`, φ² weights, KWC couple stress `kwcStress` |
| 843–878 | internal force | `gaussInt` (f_int = dE/dU at the Gauss points, energy `gaussE`), `residual` |
| 879–995 | exact Gauss Hessian | `packPd`, `gaussTangent` (one-form `vJg` assembly), `consistentTangent` |
| 996–1072 | guards and energies | stage breadcrumbs `mark` (`<out>_stage.txt`, `-trace`), NaN guards `nanE/nanV`, `refreshTangent`, `elasticEnergy`, `storedEnergy`, kinetic energies, heat |
| 1073–1139 | output | P0 → P1 projection, `writeVTK` (VTU + PVD; field names unchanged) |
| 1140–1157 | GND | `curlFpG`: G = F_p Curl F_p, ρ_G = \|G − G_R\|/b |
| 1158–1202 | time-loop setup, linear solve | step counts (`nrelax + tend/dt`), mobile-density update `mobUpdate`, `linSolve` |
| 1203–1270 | global Newton | `newtonSolve` (full Newton, line search on \|R\|, stagnation exit), plate work ledger `plateLog`, `reactionNow` |
| 1271–1369 | checkpoints and restart | atomic generations `saveCkF/saveCkA` with manifest and `.gnd/.plate/.rhoM` side files; restart reader (a relaxed checkpoint without side files is accepted) |
| 1370–1405 | run metadata | `<out>_meta.json`, CSV headers |
| 1406–1639 | time loop | per step: φ weights and drag → predictor (U + dt V, plate) → Newton (retry with fresh tangent, `_NOCONV.txt`) → Newmark commit and plate work → dissipation → φ solve (primal-dual active set) → SSD storage/recovery, GND, mobile density → heat. At `n == nrelax`: relaxed checkpoint `<out>_relaxed_ck.txt` and `_relaxed.vtu`, then the end-of-relaxation reset (V, A, T, ρ_S, θ_R, G_R). Output: log/thermo/diag rows, VTU frames, checkpoints, rolling checkpoints, TSTOP and wall-limit stops |

## Scripts (call graph; what each does and how to use it: README.md)

```
submit_ladder.sh ─┬─ run_relax.sbatch ── run_relax.sh ── run.sh (RELAXONLY=1) ── tag.sh
                  ├─ submit_speeds.sh (per angle) ─┬─ run.sbatch × N (array; afterok relax, then afterany) ── run.sh ── tag.sh
                  │                                └─ post.sbatch (array) ── tag.sh, post/animate.py, grain_post.py, tau_stats.py, package.sh
                  └─ grain_summary.sbatch ── tag.sh, post/grain_summary.py
post_all.sh ─┬─ post.sbatch (array, per angle)
             └─ grain_summary.sbatch (afterany all post arrays)
resume.sh ── recover_ck.sh, then submit_speeds.sh
```

MODE (mpi | serial, default mpi) is read by every script and passed on through `sbatch --export=ALL`.

## MPI generator `mpi/make_mpi.py`

Exact text patches on the serial file (each must match the stated number of times):
PETSc `macro_ddm` setup and the overlapping decomposition; the collective wrappers (`gmaxS`, `gsumS`, `gdotS`, `gsumEs`, `isRoot`,
`syncPen`, `syncP1`) replaced by `mpiAllReduce` versions; owned-element masks in every element sum (energies, dissipation, diagnostics);
the linear solves through PETSc (`set(A, sparams = …)`; options `-petsc` (u–θ system: GMRES + ASM/ILU), `-petscfall` (MUMPS fallback),
`-petscp1` (φ and heat), set from the environment variables `PETSC`, `PETSCFALL`, `PETSCP1` by `run.sh`); the rank-count guard (anchored on the `-j`
option line); root-only output and checkpoint writes with gathered P0/P1 arrays; `"mpi"` entry in the metadata.

## Checks `checks/`

- `make_checks.py`: builds `cosserat_kwc_check.edp` = solver + `-tancheck 1|2` (tangent against central differences, FD-step scan, KWC
  tangent alone; inserted before `newtonSolve`) + `-check 1` (driving forces and internal force against energy derivatives, appended).
- `run_tangent_check.sh`: runs it for a few loaded steps from `relax_checkpoint.txt`.
- `probe_mpi.edp`, `probe_mumps.edp`: installation probes for a FreeFem module.
