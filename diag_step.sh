#!/bin/bash
#SBATCH -J diagstep
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mem=48G
#SBATCH --hint=nomultithread
#SBATCH -t 6:00:00
#SBATCH -o logs/diagstep_%j.log
# diag_step.sh -- rerun the steps after a run's latest checkpoint in a scratch folder, with one row per Newton iteration
# (<scratch>/<tag>_newton.csv) and the timing lines, optionally with other solver options. The run folder is only read.
#   usage:  [INERTIA=1] [RUNTAG=r6] ANGLE=<deg> sbatch diag_step.sh <speed m/s> <steps> [solver options, e.g. -kwctan 3 -maxit 200]
#   (or 'bash diag_step.sh ...' inside an allocation). Options given here are appended last, so they override the run's.
#   output: diag_<tag>_<job>/ (the .out file, CSVs, VTU frames); delete it afterwards.
set -uo pipefail
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
if ! command -v FreeFem++-mpi >/dev/null; then module purge; module use "$HOME/modulefiles"; module load freefem/4.15 sqlite/3.45.0; fi
UPMS=${1:?plate speed in m/s}; NS=${2:?number of steps}; shift 2
source ./tag.sh || exit 1
D=out_$TAG; B=$TAG
[ -s $D/${B}_ck.txt ] || { echo "no checkpoint $D/${B}_ck.txt"; exit 1; }
N=$(head -1 $D/${B}_ck.txt.manifest 2>/dev/null); N=${N:-?}
NR=${NRELAX:-100}
[ "$N" != "?" ] || { echo "no $D/${B}_ck.txt.manifest: cannot tell the checkpoint step"; exit 1; }
T0=$(awk -v n=$N -v r=$NR -v d=$DT 'BEGIN{printf "%.10g", (n-r)*d}')
TE=$(awk -v t=$T0 -v k=$NS -v d=$DT 'BEGIN{printf "%.10g", t + (k+0.5)*d}')
S=diag_${B}_${SLURM_JOB_ID:-$$}; mkdir -p $S logs
# the latest generation under its plain names (hard links written with every checkpoint); no manifest, so the restart reads them
for f in $D/${B}_ck.txt $D/${B}_ck.txt.*; do
  case "$f" in *.manifest*|*.g[0-9]*|*.tmp) continue;; esac
  cp -p "$f" "$S/${f#$D/}"
done
EDP=${EDP:-cosserat_kwc_v2_mpi.edp}; NP=${NP:-${SLURM_NTASKS:-16}}
MPIRUN=${MPIRUN:-mpiexec -n $NP --bind-to none}
ARGS="-pistonup $UN -pistonut $UT -ptr ${RAMP:-0.5} -pth ${HOLD:-2.0} -tend $TE -dt $DT -vtkevery 1 -ckevery 1000000 -tstop ${TSTOP:-1400} $PHYS"
XARGS="-asm 1 -slipexact 1 -localtol ${LOCALTOL:-1e-10} -betamax ${BETAMAX:-0.785398} -verbose 1 -nlogit 1"
echo "[$(date +%T)] $TAG: steps $((N+1))..$((N+NS)) (t $T0 -> $TE ns) from $D, extra options: ${*:-none}; output $S/" | tee $S/$B.out
$MPIRUN FreeFem++-mpi -nw -v 0 $EDP $ARGS $XARGS -npexpect $NP -out $S/$B -restart $S/${B}_ck.txt "$@" >> $S/$B.out 2>&1
echo "[$(date +%T)] exit $?" | tee -a $S/$B.out
grep -h "not converged\|retry\|NOCONV\|line search failed\|^done\." $S/$B.out | head -20
echo "per-iteration rows: $S/${B}_newton.csv"
