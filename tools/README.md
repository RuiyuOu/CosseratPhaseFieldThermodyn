# tools/ (not used by the runs)

- `kwc_fit/`: fit of the KWC constants to the grain-boundary energies of Tschopp et al. (2015): `kwc_fit_tschopp.py`, `tschopp_grab.py`
  (digitised data), `kwc1d.py` (1D KWC profiles), `README_tschopp.md`. Run from this folder; outputs go to the current directory.
- `debug/`: `make_dbg.py` writes a debug build of the MPI solver (`python3 tools/debug/make_dbg.py cosserat_kwc_v2_mpi.edp
  tools/debug/cosserat_kwc_v2_dbg_mpi.edp`, from the repository root). For rerunning steps of a run use `diag_step.sh` (root).
  `make_dbg.py` was written for version 7.3; some of its anchors no longer match the 8.0 solver (it stops and names the one).
