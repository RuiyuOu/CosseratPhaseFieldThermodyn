# Stable code + dual coupling patch (2026-11-02)

Unzip on top of `hpc_dynamic/` (it only adds files and replaces scripts; existing runs and tags with default settings are unaffected).

## Files
| file | what |
|---|---|
| `cosserat2d_mob_stable.edp` | stable serial solver: previous model + `-mucdual 1` + `-kwcediff 0` (no symfold / KWC-field code) |
| `cosserat2d_mob_stable_mpi.edp` | its MPI/PETSc version (generated: `python3 mpi/make_mpi.py cosserat2d_mob_stable.edp cosserat2d_mob_stable_mpi.edp`) |
| `mpi/make_mpi.py` | generator, now for the stable and the experimental code |
| `run_piston_mob.sh`, `run_piston_mob_mpi.sh` | `CODE` (default `stable`), `MUCDUAL`, `DUALR`, `KWCEDIFF`, `PETSCDUAL` |
| `run_all_A_mpi.sh` | exports the new variables; `ARRAY` (default `0-5`) selects ladder levels |
| `run_post_piston_mob*.sbatch` | know the new tag suffixes (pass the same variables as for the run) |

## New options
* `MUCDUAL=1` (`-mucdual 1`): rotational compatibility a = skew(Ubar^e) = 0 imposed exactly by a P0 multiplier (the skew stress)
  instead of the penalty mu_c; monolithic inf-sup stable mixed Newton on [P2, P2, P1, P0] (rot u = div u^perp: Stokes P2-P0 pairing),
  exact consistent tangent incl. the local slip/spin response to the multiplier. Serial: UMFPACK; MPI: one PETSc Mat,
  `PETSCDUAL` (default `-ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps`). Tag `_du`.
  `DUALR=r` (default 0): optional augmented-Lagrangian stabilisation (GPa), tag `_dr<r>`. Multiplier side file `<checkpoint>.lam`.
  Starting from a penalty checkpoint (out_A relaxed state) the multiplier starts at 0: the first step needs a few line-search halvings.
* `KWCEDIFF=0` (`-kwcediff 0`): drop the eps^2 |grad theta|^2 term of the KWC energy (identical to `-kwce 0`). Tag `_ne`.

## Verification (local, bicrystal)
* defaults: stable == previous solver, byte-identical logs and VTK (minus the experimental fields theta_sym_deg / theta_kwc_deg)
* block tangent vs central FD: 7e-9 (momentum and constraint rows); constraint |a| ~ 1e-14; quadratic Newton convergence
* MPI 2 ranks == serial: identical tangent check and relaxation, loaded steps agree to 7 digits;
  restart MPI checkpoint (+ .lam) -> serial reproduces the uninterrupted run; penalty path MPI == serial
* `-kwcediff 0` == `-kwce 0` (identical logs)

## Update (2026-11-02, later)
* `KWCREG=1` (`-kwcreg 1`): Kobayashi & Giga (1998) regularisation of the singular term, (1/gamma) ln cosh(gamma |grad theta|),
  instead of sqrt(|grad theta|^2 + delta^2) - delta; `KWCGAM` = gamma (default 1/delta = 200). Tags `_kg`, `_kg<gamma>`.
* `KWCE` = eps of the eps^2 term (default 0.004 um, tag `_e<eps>` otherwise).
* Dual restart from a penalty checkpoint without `.lam`: the multiplier starts at the penalty skew stress 2 mu_c a.
* Checks: defaults byte-identical; KG tangent vs FD 6.7e-9; MPI 2 ranks == serial.
