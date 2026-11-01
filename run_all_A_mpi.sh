#!/bin/bash
# MPI/PETSc version of run_all_A.sh: full variant-A piston ladder (u_p = 25 50 100 150 250 400 m/s) with cosserat2d_mob_mpi.edp,
# solver-default time step DT=0.01 ns and NREF=20 (NREF is unused by full Newton), followed by the post-processing array.
#   NP ranks per level (default 8; sbatch -n is set to NP), MEM per job (default 24G); tags get the suffix _mpi (serial runs untouched).
#   physics defaults: DRAGRHO=3 VC=2 RECOV=2 CD=100 EOS=vinet VOLTAN=1 DISSIP=1 CDSTORE=1 RELAXFIX=1 ONFAIL=retry CKROLL=20
#   outputs: out_pistM3v2RcV[_<RUNTAG>]_mpi_up<speed>ms/   animations: deliver_pistM3v2RcV[_<RUNTAG>]_mpi_up<speed>ms/
# usage: bash run_all_A_mpi.sh [N_chained_arrays=10]   e.g.  RUNTAG=gm NP=8 bash run_all_A_mpi.sh
#   each array job has 24 h; a level that is complete, hit TSTOP, stalled or stopped on a non-converged step exits at once,
#   so extra links cost nothing. Override any default on the command line, e.g.  EOS=svk bash run_all_A.sh
set -euo pipefail
cd "$(dirname "$0")"
N=${1:-10}
export DT=${DT:-0.01} NREF=${NREF:-20} NP=${NP:-8}
MEM=${MEM:-24G}
export TEND RAMP HOLD VTKEVERY TSTOP GNDPHI2 RECOV CA CD DIAG DRAGRHO VC RHOMSAT C2M RHOFCUT RHOMMIN RELAXFIX \
       SPINMAX SPINBRANCH BCKIN PREDICTOR ANGLE ONFAIL MAXITRETRY CKROLL EXACTEXP DISSIP CDSTORE SOFTMIN EOS KP K0 VOLTAN NEWTON TOLRES RESFLOOR LSMAX MAXIT NLOCAL TANGENT POLISH FDK FDQ LINSOLVER LINTOL PRECEVERY TANEVERY SLIPPRED RUNTAG SHARETAG PETSC PETSCFALL PETSCP1 MPIRUN
[ -s cosserat2d_mob_mpi.edp ] || { echo "missing cosserat2d_mob_mpi.edp"; exit 1; }
[ -s out_A/A_dyn_relaxed_ck.txt ] || { echo "missing out_A/A_dyn_relaxed_ck.txt (run the variant-A relaxation first)"; exit 1; }
mkdir -p logs
prev=$(sbatch --parsable --export=ALL -n $NP --mem=$MEM run_piston_mob_mpi.sbatch)
echo "array 1/$N: $prev"
for i in $(seq 2 "$N"); do
  prev=$(sbatch --parsable --export=ALL --dependency=afterany:$prev -n $NP --mem=$MEM run_piston_mob_mpi.sbatch)
  echo "array $i/$N: $prev"
done
post=$(sbatch --parsable --export=ALL --dependency=afterany:$prev run_post_piston_mob_mpi.sbatch)
echo "post-processing array: $post (starts after $prev)"
echo "outputs: out_pistM3v2RcV[_dt${DT}_nr${NREF}][_${RUNTAG:-}]_mpi[_a${ANGLE:-45}]_up<speed>ms/   (check progress: tail -1 out_*/*_log.csv)"
