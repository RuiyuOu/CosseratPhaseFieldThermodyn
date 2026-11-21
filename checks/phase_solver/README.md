# Phase-field solver benchmarks (version 7.6)

Records of the comparison in the production notes, Section 6.13. Run from the repository root (they read `mesh.msh` and
`grain_angles_rad.txt`).
- `kkt2.py`: Python port of the phase step (matches the solver to 1.3e-11) with six solvers of the KKT system with explicit multipliers.
- `bench.edp`: FreeFEM timing of the linear solves of the active-set step (needs the standalone phase plugin it loads; mode in the last
  argument).
- `ind.py`: independence analysis of the periodic matrix check from the FEAT lines of a run with `-asmcheck 1`.
