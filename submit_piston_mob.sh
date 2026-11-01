#!/bin/bash
# Mobile-density drag piston ladder (variant A x 6 speeds; defaults DRAGRHO=3, VC=2, RECOV=2 -> tag pistM3v2R), 30 ns each, plus N-1 follow-up arrays (afterany) for the 24 h limit.
# usage: bash submit_piston_mob.sh [N_arrays=6]    optional env: TEND, RAMP, HOLD, VTKEVERY, TSTOP, DRAGRHO, VC, RHOMSAT, C2M, RHOFCUT, RHOMMIN, RELAXFIX, DT (0.01), NREF (20), SPINMAX (0=off), SPINBRANCH (pi/4), ONFAIL (retry), MAXITRETRY (80), CKROLL (20), EXACTEXP (1), DISSIP (1), CDSTORE (1), SOFTMIN (0.01), EOS (vinet), KP (5.0), K0, VOLTAN (1), NEWTON (full), TOLRES, RESFLOOR, LSMAX, MAXIT, GNDPHI2 (default 1), RECOV (default 2 = Ask recovery), CA, CD (default 100 with RECOV=2, else 10), DIAG
# needs out_A/A_dyn_relaxed_ck.txt (same relaxed start as the original A ladder)
N=${1:-6}; mkdir -p logs; export TEND RAMP HOLD VTKEVERY TSTOP GNDPHI2 RECOV CA CD DIAG DRAGRHO VC RHOMSAT C2M RHOFCUT RHOMMIN RELAXFIX DT NREF SPINMAX SPINBRANCH BCKIN PREDICTOR ANGLE ONFAIL MAXITRETRY CKROLL EXACTEXP DISSIP CDSTORE SOFTMIN EOS KP K0 VOLTAN NEWTON TOLRES RESFLOOR LSMAX MAXIT NLOCAL TANGENT POLISH FDK FDQ LINSOLVER LINTOL PRECEVERY TANEVERY SLIPPRED RUNTAG
[ -s out_A/A_dyn_relaxed_ck.txt ] || { echo "missing out_A/A_dyn_relaxed_ck.txt"; exit 1; }
prev=$(sbatch --parsable --export=ALL run_piston_mob.sbatch)
for i in $(seq 2 $N); do prev=$(sbatch --parsable --export=ALL --dependency=afterany:$prev run_piston_mob.sbatch); done
echo "mobile-drag piston ladder submitted; last array job $prev"
echo "outputs: out_pistM<DRAGRHO>v<VC>[G][R]_up<m/s>ms/ (default pistM3v2RcV)  (frames every 0.1 ns, <tag>_diag.csv = worst element per step); a STALLED file marks a failed level"
