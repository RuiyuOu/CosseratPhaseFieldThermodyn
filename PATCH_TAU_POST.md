# Patch 2026-11-09: everything after the 3-point (-nq 3) quadrature work

Unzip on top of the project root (the folder that contains `hpc_dynamic/` and `docs/`). Existing runs and default solver behaviour are unchanged.

| file | what |
|---|---|
| `hpc_dynamic/cosserat2d_mob_stable.edp`, `cosserat2d_mob_stable_mpi.edp` | solver state after the quadrature work (memory-split Gauss tangent). NOTE: the MPI file was generated before the last serial edit and was not regenerated; its tangent assembly is structured differently (`gaussTangentPoint`) from the serial one (`gtBlocks`/`gtAccum`). Regenerate with `python3 mpi/make_mpi.py cosserat2d_mob_stable.edp cosserat2d_mob_stable_mpi.edp` if you want them identical; not tested on Longleaf. |
| `hpc_dynamic/quadtest/` | quadrature comparison harness (`run_quad.sh`, `chain.sh`, `compare_quad.py`, `args_base.txt`, bicrystal mesh) |
| `hpc_dynamic/post/tau_stats.py` | tau_1..3, tau_c, T, pressure, heating-jump and mismatch-beta statistics over the whole run |
| `hpc_dynamic/post/package.sh` | optional 6th argument: tau_stats output folder -> `deliver_<tag>/tau/` |
| `hpc_dynamic/run_post_piston_mob.sbatch`, `..._mpi.sbatch` | call tau_stats.py after the grain post-processing (env `TAUWIN0 TAUWIN1 TAUMARK TAUEVERY`) |
| `hpc_dynamic/README.md`, `CODE_MAP.md` | documentation of the above |
| `docs/quadrature_test.{md,pdf}` | quadrature test description |
| `docs/cosserat_kwc_production.{md,pdf}` | production-configuration implementation document (incl. sections 3.5-3.7) |
| `docs/cosserat_kwc_theory_3d.{md,pdf}` | general 3D theory document |

Testing status: `tau_stats.py` and `package.sh` tested locally on synthetic VTU series built from one real snapshot; the sbatch files are syntax-checked only; nothing here has been run on Longleaf.
Usage: `python3 post/tau_stats.py out_<tag> --window 0.5 30 --mark 6 --log out_<tag>/<tag>_log.csv`
