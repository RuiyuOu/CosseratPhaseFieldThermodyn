# AD tangent for the 2D solver (default)

The 2D solver runs by default (`TANGENT` unset or `ad`) with the element-map derivatives (consistent tangent, local slip/spin Newton) computed exactly by forward-mode
automatic differentiation in the plugin `kwc2d.cpp` (`dual.hpp`), instead of central differences (`fdMap`). Nothing else changes: the AD file
`cosserat_kwc_ad.edp` is generated from `cosserat_kwc.edp` by `make_ad.py` (and its MPI version by `mpi/make_mpi.py`), so the stable solver stays
the single source. Default: `ad`. If the plugin cannot be built (no `ff-c++`, compiler error) a run with `TANGENT` unset falls back to `fd` with a
WARNING in the run's `.out`; `TANGENT=ad` fails instead, `TANGENT=fd` forces the original finite-difference solver. Results agree, so a run
may be resumed with either tangent.

| file | what |
|---|---|
| build_ad.sh | builds kwc2d.so (into the repo root, where FreeFem loads it) and generates cosserat_kwc_ad.edp / cosserat_kwc_ad_mpi.edp when stale; run.sh calls it |
| make_ad.py | cosserat_kwc.edp -> cosserat_kwc_ad.edp (fdMap -> AD; `-adcheck 1` = AD vs central differences on every element, then stop) |
| kwc2d.cpp, dual.hpp | the plugin (element map of trialFp / finiteStress / eosStress / taus) and the dual-number type |
| ad_compare.sbatch | fd vs ad on 16 MPI ranks from relax_checkpoint.txt (RUNTAG fdcmp / adcmp), then compare_runs.py |
| compare_runs.py | per-step differences of two runs (log columns, iterations, last common VTU, wall times) |

```
sbatch ad/ad_compare.sbatch                         # UPMS=400 ANGLE=45 TEND=1.0 by default; result in logs/adcompare_<job>.log
MODE=mpi bash run.sh 400                            # default: AD tangent
TANGENT=fd MODE=mpi bash run.sh 400                 # the finite-difference tangent
```
Local check (serial, 2 relaxation + 3 loaded steps at 400 m/s): identical iteration counts, u within 1.4e-14 um, phi and T identical,
wall time 446 s -> 365 s. -adcheck: AD vs central differences 1e-7 relative (the FD error), values 3e-13.

Cluster check (Longleaf, 16 MPI ranks, 400 m/s at 45 deg, 100 loaded steps from the relax checkpoint, 2026-11-11): 607 vs 607 Newton iterations,
0 steps differ, every logged column identical (E_elastic, KE, Q, slip, dT, dtheta, rhoG, phi, audit); wall 1469 s (fd) -> 1206 s (ad), -18%.
-adcheck on 58391 elements: AD vs central differences 2e-5 of 250 (the FD error), plugin vs solver values 2.8e-13, Jacobians 0.65 s (FD) vs 0.12 s (AD).
