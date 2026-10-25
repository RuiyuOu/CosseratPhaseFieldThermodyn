# HPC bundle: dynamic Cosserat-KWC 15 ns, variants A and F
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
Post-process one level: cd post; BASE=lad_P1000MPa_S100MPa NP=4 python3 animate.py ../out_lad_P1000MPa_S100MPa ../anim_lad_P1000MPa_S100MPa
