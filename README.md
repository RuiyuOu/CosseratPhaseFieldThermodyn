# Cosserat–KWC plate impact (FreeFem++)

Dynamic finite-deformation plane-strain Cosserat crystal plasticity with the KWC orientation phase field, on an 8-grain
[1-10] Rice planar-triple-slip polycrystal loaded by a moving plate. Units: µm, ns, GPa, K.

## Production model

- **Kinematics:** F = R(θ) Ū F_p, F_p updated by the exact exponential of the slip (3 systems) and grain-boundary spin increments.
- **Elasticity:** Vinet volumetric law, isochoric cubic Green-strain energy and Cosserat penalty 2μ_c a², integrated at the 3 edge midpoints of each P2 triangle.
- **KWC:** ψ = f₀[(1−φ)²/2 + a²|∇φ|²/2 + s g(φ) N(∇θ)], with g = −2[φ + ln(1−φ)] and N = √(|∇θ|² + δ²) − δ.
- **Plasticity:**
  - thermal activation in series with relativistic drag;
  - mobile density ρ_M evolved, with drag density min(ρ_M, ρ_S + ρ_G);
  - GND density from G = F_p Curl F_p;
  - sweeping recovery (Ask et al. 2019), with storage capped by the dissipation.
- **Loading:** left edge clamped; the right edge is the plate (ramp 0.5 ns, hold 2 ns); top and bottom free. A plate speed V at ANGLE splits into compression V cos(ANGLE) along −x and shear V sin(ANGLE) along +y.
- **Numerics:**
  - Newmark average acceleration;
  - full Newton with the exact Gauss-point tangent and a line search;
  - backward-Euler local slip and spin update;
  - staggered phase-field solve (primal-dual active set) and heat solve.
- **Start:** 100 zero-load relaxation steps (φ starts at 0.999), then the load. DT 0.01 ns, t_end 30 ns.

Parameter values are the defaults in the options table below; their derivation is in `cosserat_kwc_production.md`.

## Running on Longleaf

```bash
bash setup_py.sh                                   # once, on a login node: post-processing venv ($HOME/dynpost_venv)
bash submit_ladder.sh                              # angles 0 45 90 x speeds 25 50 100 150 250 400 m/s (16 ranks, 48 GB per run)
ANGLE=0 bash submit_speeds.sh                      # one angle, all speeds
ARRAY=5 ANGLE=90 bash submit_speeds.sh             # one angle, one speed (5 = 400 m/s)
RUNTAG=s12 KWCS=1.2 bash submit_ladder.sh          # changed parameters need a new RUNTAG
ANGLE=45 bash resume.sh                            # resubmit after crashes
squeue -u $USER;  tail -qn 1 out_a*_up*ms/*_log.csv # queue and progress (column 2 = time in ns)
```

- Every script runs the MPI solver unless `MODE=serial` is set, for example `MODE=serial ANGLE=0 bash submit_speeds.sh`. Serial runs use 1 core, DT 0.005 ns and the tag suffix `_serial`.
- The submit scripts queue the relaxation job first. It makes `relax_checkpoint.txt` if it is missing and does nothing if the file exists with the same parameters; every run starts from it.
- Rerunning a submit script resumes unfinished runs and skips finished ones. Each 24 h array link continues from the last checkpoint, and links of finished runs exit at once.
- `cosserat_kwc_mpi.edp` is generated from the serial solver: `python3 mpi/make_mpi.py cosserat_kwc.edp cosserat_kwc_mpi.edp`. Never edit it by hand; the submit scripts regenerate it if it is missing.

## Scripts

**Submitting (login node):**

| Script | What it does | Usage |
|---|---|---|
| `submit_ladder.sh` | The full study: one relaxation job, then `submit_speeds.sh` for each angle, then `grain_summary.sbatch` after all post-processing. | `bash submit_ladder.sh [N]`; env `ANGLES` ("0 45 90"), `SERIAL_ANGLES=1` (one angle at a time, at most 96 cores) |
| `submit_speeds.sh` | One angle, all speeds: a relaxation job, N chained 24 h arrays of `run.sbatch` (N = 10, serial 20), then `post.sbatch`. | `ANGLE=0 bash submit_speeds.sh [N]`; env `ARRAY` (0–5 = 25…400 m/s), `NP` (16), `MEM` (48G) |
| `post_all.sh` | Post-processing of every angle now (from the frames written so far), then the grain summary once those jobs finish. | `bash post_all.sh`; env `ANGLES` ("0 45 90"), `ARRAY`, `FORCE=1` (redraw all animation frames), `RUNTAG`/`DT`/`MODE` as for the runs |
| `resume.sh` | After crashes: clears the stop markers of one angle's runs, rolls back a TSTOP above 5000 K (blow-up) with `recover_ck.sh`, then resubmits with `submit_speeds.sh`. Finished runs and genuine TSTOPs are left alone. | `ANGLE=0 bash resume.sh`; env `N` (20), `TBLOW` (5000), `FORCE=1` (even if a run looks active) |
| `recover_ck.sh` | Rolls one run back to its newest rolling checkpoint before the bad step and clears its stop markers. | `bash recover_ck.sh out_<tag> [...]` |
| `setup_py.sh` | One-time: creates `$HOME/dynpost_venv` with numpy, scipy, matplotlib, meshio, pillow, pandas. | `bash setup_py.sh` |

**Slurm jobs** (queued by the submit scripts; resources set on submission):

| Job | What it does |
|---|---|
| `run_relax.sbatch` → `run_relax.sh` | Relaxation stage: makes `relax_checkpoint.txt` (output in `out_relax/`), or exits at once if it exists with the same parameters. A file made with other parameters is refused; delete it (and its `.args`) or set `RELAXCHECK=0`. |
| `run.sbatch` → `run.sh` | One speed per array task (task 0–5 = 25 50 100 150 250 400 m/s), 24 h per link. |
| `post.sbatch` | One speed per array task: animation frames (`anim_<tag>/`) and the delivery folder `deliver_<tag>/` (gifs, contact sheets, png zips, data zip, `grains/` grain history, `tau/` stress statistics). Run it any time to post-process the frames written so far: `ANGLE=0 sbatch --export=ALL --array=0-5 post.sbatch`. Env: `TAUWIN0 TAUWIN1` (window, 0.5–30 ns), `TAUEVERY`, `FORCE=1` (redo frames). |
| `grain_summary.sbatch` | Grain refinement and misorientation summary over angles × speeds → `grain_summary[<suffix>]/`. Env: `ANGLES`, `UPS`, `MISOR`, `AREA`, `SMOOTH`, `TAVG`. |

**Single runs and helpers:**

| Script | What it does | Usage |
|---|---|---|
| `run.sh` | One restart-safe run, run directly or from `run.sbatch`. It checks the settings against `<tag>.args`, starts from `relax_checkpoint.txt` or resumes from the last checkpoint, picks the MPI launcher, and reports COMPLETE / INCOMPLETE / STOPPED / STALLED. | `NP=4 bash run.sh 400 0.5` (speed m/s, t_end ns [, ramp, hold, vtkevery]); `MODE=serial bash run.sh 400 0.5` |
| `run_relax.sh` | The relaxation stage on its own (as in `run_relax.sbatch`). | `NP=4 bash run_relax.sh`, `MODE=serial bash run_relax.sh` |
| `tag.sh` | Sourced by every script, never run: the run tag and the full option list `PHYS`. | — |
| `checks/run_tangent_check.sh` | Serial self-test of the tangent and energy derivatives (see Checks). | `bash checks/run_tangent_check.sh [steps]` |

`MODE`, `ANGLE`, `RUNTAG`, `DT` and the physics variables must be the same for the runs, the post-processing and `resume.sh`. They are, as long as everything is submitted with the same environment.

## Outputs

Every column and field of these files, with its unit and formula, is described in `DATA_FIELDS.md`.

| Path | Contents |
|---|---|
| `out_<tag>/` | `<tag>_log.csv`, `_thermo.csv`, `_diag.csv`, VTU frames + `<tag>.pvd`, `_meta.json`, checkpoints (`_ck.txt` + manifest, `_ckroll{0,1}.txt`), `<tag>.args` |
| `deliver_<tag>/` | results: per-field gif, contact sheet and png zip, `<tag>_data.zip` (raw run data), `grains/` (grain count, size and misorientation history), `tau/` (resolved-stress, temperature and pressure statistics) |
| `anim_<tag>/` | animation frames (reused by later post runs; `FORCE=1` redoes them) |
| `grain_summary[<suffix>]/` | summary over angles and speeds |

- **Tag format:** `a<ANGLE>_up<speed>ms[_dt<DT>][_<RUNTAG>][_serial]`.
- **Resume guard:** a resume with other settings than those in `<tag>.args` is refused, and the message names the option.
- **Stop markers:** `_TSTOP.txt` (T above 1400 K), `_NOCONV.txt`, `_NAN.txt`, `STALLED`. `resume.sh` clears them, and rolls back a TSTOP above 5000 K as a numerical blow-up.

## Options

`tag.sh` passes all of these to the solver. To change one, set the environment variable of the same name in capitals (for example `KWCS=1.2`).

| Group | Option (default) |
|---|---|
| KWC | `-f0` 0.005944, `-kwca` 0.10, `-kwcs` 1.05, `-kwcd` 0.005 |
| phase-field kinetics | `-etaphi` 0.1, `-qphi` 100e3, `-etaphimin` 1e-2 |
| GB spin, recovery | `-etastar` 1, `-cd` 100, `-ca` 3.16227766 |
| dislocations | `-vc` 2, `-rhomsat` 100, `-c2m` 10, `-rhofcut` 1e4, `-rhommin` 1 |
| Cosserat, elasticity | `-muc` 75.4, `-j` 0.01, `-kp` 5.0, `-k0` 0, `-softmin` 0.01 |
| relaxation | `-nrelax` 100 |
| Newton | `-maxit` 40, `-maxitretry` 80, `-tolres` 1e-6, `-resfloor` 1e-8, `-lsmax` 6, `-stagres` 100, `-jmin` 0.1, `-tanevery` 1, `-precevery` 0, `-lintol` 1e-10 |
| local update | `-nlocal` 30, `-polish` 20, `-fdk` 1e-6, `-fdq` 1e-8, `-spinbranch` 0.785398 |
| phase-field solve | `-phitol` 1e-10, `-phikkt` 1e-8, `-phimaxit` 40 |
| output | `-diag` 1, `-ckroll` 20 |

Run settings (environment): `ANGLE` (45), `TEND` (30), `RAMP` (0.5), `HOLD` (2.0), `DT` (0.01), `VTKEVERY` (10), `TSTOP` (1400), `RUNTAG`, `NP`, `MEM`, `ARRAY`.

## Mesh and grain angles

`mesh.msh` and `grain_angles_*.txt` are made from a Neper tessellation, whose files are kept in `mesh/`:

```bash
# 1. tessellation (Neper, built from github.com/neperfepx/neper): 8 grains in a 10 x 10 um square, seed 1
neper -T -n 8 -dim 2 -id 1 -reg 1 -domain "square(10,10)" -oridescriptor euler-bunge -format tess,geo,ori -o MPTD_2D
# 2. FreeFem mesh and grain angles
python3 mesh/make_mesh.py mesh/MPTD_2D.geo mesh/MPTD_2D.ori --out mesh_build
cd mesh_build && FreeFem++ -nw make_mesh.edp && cd ..
cp mesh_build/mesh.msh mesh_build/grain_angles_deg.txt mesh_build/grain_angles_rad.txt .
```

- **Mesh:** FreeFem `buildmesh` with one straight border per line of `MPTD_2D.geo`. Grain boundaries get segments of about 0.05 µm and boundary label 100 + line number; the outer sides get about 0.1 µm segments and labels 1 (left), 2 (bottom), 3 (right), 4 (top). Each triangle's region is its grain number, the `.geo` surface number, which is also the row of `MPTD_2D.ori`.
- **Grain angles:** θ = Euler–Bunge φ₁ reduced modulo 180° into [−180°, 0). With [1̄10] out of plane the lattice repeats every 180°.
- **Check:** these steps reproduce the repository's `mesh.msh` (29,399 vertices, 58,391 triangles, 1,227 boundary edges; same triangles and labels, vertices within 5×10⁻⁷ µm) and its angle files exactly.

## Checks

`bash checks/run_tangent_check.sh` checks the tangent and the energy derivatives against finite differences. `checks/probe_mpi.edp` and `checks/probe_mumps.edp` test a FreeFem module's MPI/PETSc and MUMPS.
