# Post-processing changelog

Post-processing is versioned separately from the solver, which is frozen at version 8.0 (`SOLVER_SHA256`). Nothing in `post/` or in
`post.sbatch`, `post_all.sh` and `grain_summary.sbatch` can change a solution; record every change of a figure, a statistic or a default
here, with the date, so that a figure in the manuscript can be traced to the version that drew it.

| Date | Change |
|---|---|
| 2026-11-20 | Baseline at the solver freeze: `animate.py`, `grain_post.py` (grain recovery, misorientation; `MINAREA`), `grain_count.py`, `grain_summary.py`, `tau_stats.py`, `wavepost.py`, `sliplaw_compare.py`, `compare_fields.py`, `package.sh`. |
| 2026-11-20 | drx_evidence.py: DRX tests per run (defects and GOS of new vs remnant grains, Derby, similitude, Hughes scaling, HAGB fraction, nucleation sites, stored fraction, GB-spin share); output deliver_<tag>/drx/ |
| 2026-11-20 | drx_evidence.py: DRX tests per run (defects and GOS of new vs remnant grains, Derby, similitude, Hughes scaling, HAGB fraction, nucleation sites, stored fraction, GB-spin share); output deliver_<tag>/drx/ |
