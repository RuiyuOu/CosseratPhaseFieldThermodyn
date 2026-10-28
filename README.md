The deliver_pist_up250ms directory has the gif animations and graphics used in the write up. the anim directories have pngs of each snapshot. 
the out_ directories have the data for each snapshot. Please do not git clone the entire thing, there are GBs of data files here. 

# HPC bundle: dynamic Cosserat-KWC
Contents: cosserat2d.edp (finite-deformation solver), mesh_A.msh, mesh_F.msh, run_dyn.sh, run_dyn_{A,F}.sbatch, run_post.sbatch,
submit_all.sh, post/ (wavepost.py, animate.py, package.sh, SMOOTHING.md).
Use: unzip, cd into the folder, `bash submit_all.sh [N]` submits a chain of N (default 4) restart-safe jobs per variant: each resumes from the last checkpoint, so no manual resubmission is needed. Manual: `sbatch run_dyn_A.sbatch` again resumes too.
Single core per run (the solver is serial FreeFem++; ~11-25 s/step locally => ~5-10 h; HPC cores should be faster). Peak RAM ~3 GB (8G requested).
Output: out_A/A_dyn_s{N}.vtu and out_F/F_dyn_s{N}.vtu every 5 steps = 0.05 ns (~300 frames each, all fields), .pvd, _log.csv, _meta.json, _ck.txt.
Restart-safe: rerunning a sbatch resumes from the last checkpoint (every 50 steps = 0.5 ns).
Post: run_post.sbatch makes 14 animations x {GIF, PNG zip, contact sheet} per variant plus raw-data zips in deliver_A, deliver_F.
Edit the module lines (freefem/4.15 sqlite/3.45.0) / python module name to match the cluster. Solver loads plugin iovtk (FF_LOADPATH is set if found).
Solver fixed vs older bundle: string literals split to compile on 4.13/4.15.

## Post-processing on the HPC (after the runs finish)
1. One-time: `bash setup_py.sh` on a login node (creates ~/dynpost_venv with numpy, scipy, matplotlib, meshio, pillow).
2. `sbatch run_post.sbatch` (or it is queued automatically by submit_all.sh). Needs out_A/ and/or out_F/ (the *_s*.vtu frames + .pvd) in the bundle folder.
3. Results: anim_{A,F}/ (14 animations x GIF and png/), deliver_{A,F}/ (GIFs, contact sheets, PNG zips, data zip, SMOOTHING.md).
Single variant / test: cd post; BASE=A_dyn NP=4 python3 animate.py ../out_A ../anim_A
Cost: ~300 frames x 14 plots; a few minutes per variant on 4 cores. Memory ~2-4 GB per worker.

## Compression + shear test (default: 10 GPa compression + 1 GPa shear, variant A, resumes from out_A/A_dyn_relaxed_ck.txt)
`bash submit_comp.sh` (chain of 6 restart-safe jobs). Env overrides: AMP, TEND, TSCALE, TAG. Output (auto tags): out_cs10_1/cs10_1_*.vtu for 10 GPa compression + 1 GPa shear, out_comp10/comp10_*.vtu for SHEAR=0 (pure compression); the original shear run stays in out_A/A_dyn_* (every 0.05 ns).
Load: traction (-AMP, +SHEAR) GPa on the right edge (label 3) with the left edge clamped (new solver option -dir combo -shear; SHEAR=0 gives pure compression); envelope = 2 ns smooth ramp, hold to 5 ns, 2 ns unload (TSCALE stretches/compresses it).
Watch min_phi / max_rhoG / mean_rhoS in the *_log.csv in the case folder for nucleation (see notes).

## Load ladder (compression P + shear 0.1 P; all levels in parallel)
`bash submit_ladder.sh [N]` -> 7 array tasks: P = 100, 200, 500, 1000, 2000, 5000, 10000 MPa (shear 10 ... 1000 MPa), each from out_A/A_dyn_relaxed_ck.txt.
Outputs out_lad_P<MPa>MPa_S<MPa>MPa/. Chain of N arrays gives restart after 12 h limits. A level that blows up (NaN) writes out_<tag>/STALLED with the end of the solver log and stops; other levels are unaffected.
Compare levels with min_phi, max_rhoG_um2, mean_rhoS_um2 in each log (nucleation shows as phi dropping where rhoS is high; threshold rho ~ 3.8e3 um^-2).
Post-process all levels: `sbatch run_post_ladder.sbatch` (array, one task per level; outputs deliver_<tag>/ with GIFs, contact sheets, PNG zips, data zip; needs setup_py.sh once).
Post-process one level: cd post; BASE=lad_P1000MPa_S100MPa NP=4 python3 animate.py ../out_lad_P1000MPa_S100MPa ../anim_lad_P1000MPa_S100MPa

## Nucleation probe
`cd post; python3 nucleus_probe.py ../out_lad_P10000MPa_S1000MPa lad_P10000MPa_S1000MPa ../nuc_P10000 [dphi_thr=0.7] [dtheta_deg=2] [frame_step=1]`
-> nuc_P10000_nucleus.csv (per frame: area of new phi-drop regions where phi was ~1 at the reference frame, area of lattice-orientation pockets
(misrotation vs grain median developed since the first frame), largest region size in nm, position, mean T and rhoS) and nuc_P10000_nucleus_last.png.
Baseline: the 15 ns shear run gives zero new phi-drop and zero pocket area (negative control).

## Confined piston ladder (cold-spray analogue), 30 ns, variants A and F  -- NEW solver option -dir piston
Set-up: wall u_x=0 (left), rollers u_y=0 (top and bottom, exactly 0), right edge prescribed u_x = -d(t); d(t) = integral of u_p*e(t), e = smooth quintic ramp up (-ptr, default 0.5 ns), hold (-pth, default 2 ns), ramp down;
piston travel = u_p*(ptr+pth) (u_p=250 m/s -> 0.625 um = 6.25 %); afterwards the piston stays put and the block evolves to 30 ns. Arguments: -pistonup (um/ns = 1000 m/s), -ptr, -pth, -tstop (K).
`bash submit_piston.sh [N]`: 12 array tasks = variants A (tasks 0-5, out_pist_up<m/s>ms) and F (tasks 6-11, out_pistF_up<m/s>ms, no relaxation stage, starts fresh), u_p = 25, 50, 100, 150, 250, 400 m/s (impact 50 ... 800 m/s; Cu shock-fit pressure 0.9 ... 16 GPa), TEND=30 ns, dt=0.01 ns (3000 loaded steps),
frames every 0.1 ns (VTKEVERY=10 steps; ~300 frames, ~13 MB each), 24 h per job, N chained arrays (default 6) resume from checkpoints. Outputs out_pist_up<m/s>ms/.
Safeguard: -tstop 1400 K (override with TSTOP=...): a level that exceeds it (or gives NaN) writes a final frame (..._tstop<step>.vtu), a checkpoint and <tag>_TSTOP.txt and ends cleanly (counted as finished).
Notes: W_ext / audit_residual columns are not meaningful in piston mode (no traction work is accumulated). The relaxed state is reused (left wall u_y freed, rollers set u_y=0 at the first piston step).
Post: `sbatch run_post_piston.sbatch` (needs setup_py.sh once); nucleus_probe.py works on these folders too (BASE=pist_up250ms).

Variant F only: `bash submit_piston_F.sh [N]` (6 array tasks, out_pistF_up<m/s>ms/) and `sbatch run_post_piston_F.sbatch` for its post-processing.
