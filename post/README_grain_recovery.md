# Grain formation, recovery and retention diagnostics

`post/postprocess_grain_recovery.py` reads the saved VTU snapshots of one run and produces a synchronised three-panel animation, material-point and patch histories, an SSD budget per patch and a diagnostic summary. It does not rerun or change the solver, and it never writes into the run folder.

## Run

```bash
python3 post/postprocess_grain_recovery.py --config post/grain_recovery.toml \
        --run out_a90_up250ms_v2_r6 --out grain_recovery_a90_up250ms --np 8
```

- `--run`, `--base` and `--out` override the configuration; `--np` sets the number of processes; `--no-mp4` skips the animation.
- To repeat an analysis exactly, pass the `resolved_config.json` it wrote: `--config grain_recovery_a90_up250ms/resolved_config.json`. That file holds the selected points as manual points.
- TOML needs Python ≥ 3.11. On an older Python, write the configuration as JSON (same keys) or start from `resolved_config.json`.

**Dependencies:** numpy, scipy, matplotlib, pillow, and ffmpeg for the MP4. It reuses `post/grain_count.py`, which must sit in the same folder. No meshio or VTK is needed.

**Time:** about 0.3 s per production-mesh frame for the analysis and about 1 s per frame for rendering, per process. For 301 frames on 8 processes, expect a few minutes.

## What it does

1. **Frames.** It reads `<base>_relaxed.vtu` (t = 0) and every frame listed in `<base>.pvd`. Missing files, duplicate times, and gaps against the nominal spacing are reported, not filled.
2. **Material correspondence.**
   - **The mesh is a fixed reference mesh.** Every frame must have the same points and triangles as the first one; this is checked by hash, and the script stops otherwise. MPI runs write the global mesh from rank 0, so node and element indices are the same material in every frame.
   - **Points** are fixed (triangle, barycentric) positions in reference coordinates. P1 fields (θ_t, φ) are interpolated linearly; P0 fields (ρ_S, ρ_G) are taken from the containing element.
   - **Patches** are fixed sets of triangles. Their means and standard deviations are weighted by reference area. The standard deviation is spatial variability, not a confidence interval.
3. **Segmentation** uses `grain_count.segment` on the stored `theta_lattice_deg`, with the thresholds of `grain_post.py`: 2°/µm, 0.1 µm², 2° merge.
   - No wrapping and no symmetry reduction are applied to θ_t. The 90° period appears only in the misorientations grain_count reports.
   - ρ_G is not used for the segmentation.
4. **Selection** is made on the final frame and then fixed in reference coordinates.
   - **Default, mode `centers`:** one point at the centre of every final grain larger than 0.1 µm², at most 60 grains, largest first. The centre is the deepest point of the grain, the one farthest from its grain boundaries and from the domain edge. With `center_method = "centroid"` it is the area-weighted centroid instead, falling back to the deepest point when the centroid lies outside the grain.
   - **Labels:** each point is labelled `center-new` or `center-kept`. A grain counts as new if its material has rotated by at least `new_grain_dtheta_deg` since the first frame, or differs by that much from the remnant of its initial grain. The label never changes which points are tracked.
   - **Mode `auto`:** instead tracks the largest new grains, each with an interior and a near-boundary point, plus one unrefined grain.
   - **Mode `manual`:** tracks only the points given in `[[selection.points]]`.
   - **Recorded for each point:** coordinates, patch radius and area, depth and the rationale, in `selection.csv` and `summary.txt`. Final grain ids are labels, not tracked lineages.
5. **Histories** of θ_t, ρ_S, ρ_G, the forest density ρ_F = ρ_S + ρ_G (absolute and as a percentage of its own maximum over all analysed frames) and the accumulated slip rate, for every point and patch (φ is in the CSVs, not plotted), plus the distance to the instantaneous segmented boundary (a crossing explains sweeping recovery).
   - **Accumulated slip rate:** the change of the saved accumulated slip `slip_path_eff` (Σ over systems of |Δγ|, the field of the `slip_acc` animation) between consecutive snapshots, divided by their time difference, shown at the later snapshot. It is an average over the snapshot interval, not the instantaneous rate, and has no value at the first frame. A decrease (a reset of the slip path) is shown as 0 and counted in `summary.txt`.
6. **SSD budget per patch.**
   - **Law:** storage K/b √ρ dΓ, dynamic recovery 2d/b ρ dΓ, and sweeping recovery C_D ⟨dφ⟩₊ tanh(C_A²|∇θ_t|²) ρ for `-recov 2`.
   - **Where the parameters come from:** C_D, C_A and the recovery mode are read from the run's `_meta.json`, `.args` or `.out` option dump. K, d, b and the floor come from the configuration, flagged as documented constants.
   - **Estimates:** the 7.x solver saves no rate fields, so the terms are snapshot estimates. Storage and dynamic recovery use the exact cumulative `slip_path_phys`. Sweeping recovery uses the snapshot change of φ; it is a lower bound, because the positive part is taken after differencing.
   - **Comparison:** the cumulative estimates are compared with the observed change. The storage cap and the floor cannot be reconstructed.
7. **Relaxation diagnostics** from `<base>_log.csv`: the kinetic energy and Q_heat/W_plate, with their values at 13 and 30 ns. A stress measure is taken from the snapshots: the area-RMS and 99th percentile of |σ|_F.

## Outputs

| File | Content |
|---|---|
| `grain_recovery.mp4` | θ_t (twilight, fixed range), log ρ_S, log ρ_G (fixed ranges, display floor) and log accumulated slip rate (turbo, fixed range) maps in reference coordinates, with the same boundaries and markers; histories below (θ_t, ρ_S, ρ_G, ρ_F, ρ_F in % of its maximum, slip rate) with the time cursor and dashed lines at 3 and 13 ns. One saved snapshot per video frame. |
| `grain_recovery_deformed.mp4` | `deformed = true`: the same maps in the deformed geometry (x = X + u); markers and boundaries are mapped with the nodal displacement. |
| `snapshots.pdf/.png` | the same layout at 0, 3, 13 and 30 ns (where available) and at the formation frame of each tracked new grain |
| `histories.pdf/.png` | all points (solid) and patches (dotted, with the s.d. band) |
| `recovery_<P>.pdf/.png` | cumulative storage, dynamic recovery, sweeping recovery, predicted and observed Δρ_S of each patch |
| `relaxation.pdf/.png` | kinetic energy, Q_heat/W_plate, stress measure |
| `point_history.csv` | t, point id, X, Y, θ_t (raw stored), ρ_S, ρ_G, ρ_F, ρ_F in % of its maximum, φ, slip path, distance to the instantaneous boundary, instantaneous grain label and its mean θ |
| `patch_history.csv` | t, patch id, fixed reference area, weighted means and s.d., cumulative rate estimates, and their source (saved rate fields or snapshot estimate) |
| `selection.csv` | the points and patches with the rationale |
| `summary.txt` | missing snapshots, mapping assumptions, parameters and their sources, colour ranges and clipping, display and solver floors, negative densities, resets, budget mismatch, observations (numbers only) |
| `resolved_config.json` | everything needed to repeat the analysis |

## Reading the results

- **ρ_S at fixed material:** a decrease is net SSD depletion in the model. Recovery terms can be nonzero while ρ_S still rises, if storage is larger.
- **Boundary/interior structure:** low ρ_G and uniform θ_t inside a grain, with high-ρ_G bands at its edges, support a boundary/interior structure. They do not establish recrystallisation. Agreement of the GND bands with the orientation boundaries is internal consistency of coupled model fields.
- **Static recovery:** with the stated SSD law there is no static recovery. When slip and φ evolution stop, ρ_S stops changing, so continued recovery after 13 ns is not expected.
- **Stored energy:** ψ_s = ½χμb² φ ρ_S also changes with φ, so a drop in stored energy is not by itself SSD annihilation.
- **Resolution and duration:** events faster than the snapshot spacing are unresolved. Retention to 30 ns shows persistence over the simulated interval only.
- **The saved ρ_G is the measure of the code version that wrote it:** |F_p c − G_R|/b before 7.5, and the Cermelli–Gurtin |c/J_p|/b from 7.5. `summary.txt` records the version where the run logged it.

## Tests made

- **Production mesh:** 14 snapshots (0.1–3.9 ns) of the October 50 m/s run. Frame discovery reported the 27 missing files correctly, and the hash check passed. For the comparison patch, the snapshot-estimated budget gave +36.6 µm⁻² against +36.55 observed. That run formed no new grains, which the script reports.
- **Synthetic run** (coarse mesh, 31 frames, 0–30 ns): a disc rotates 10° into a new grain between 1 and 3 ns, and its interior ρ_S then falls from 60 to 20 µm⁻². The script found and tracked the disc, placed its interior, near-boundary and unrefined points, and drew the 3 and 13 ns markers. It flagged the imposed decay as unexplained by the SSD law (mismatch −9.7 µm⁻²), as it should. A rerun from `resolved_config.json` reproduced the histories exactly, and the run folders were not modified.
