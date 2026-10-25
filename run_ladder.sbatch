#!/bin/bash
#SBATCH -J ladder
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=8G
#SBATCH -t 12:00:00
#SBATCH --array=0-6
#SBATCH -o logs/ladder_%A_%a.log
# Load ladder: compression P with shear = 0.1*P, all 7 levels run at the same time (one array task each).
#   task: 0     1     2     3      4      5      6
#   P:    100   200   500   1000   2000   5000   10000 MPa     (shear 10 ... 1000 MPa)
# Each task resumes from out_A/A_dyn_relaxed_ck.txt and is restart-safe; use submit_ladder.sh for automatic follow-up arrays.
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
module purge
module use "$HOME/modulefiles"
module load freefem/4.15 sqlite/3.45.0
mkdir -p logs
LEVELS=(100 200 500 1000 2000 5000 10000)      # compression, MPa
P=${LEVELS[${SLURM_ARRAY_TASK_ID:-0}]}
S=$(( P / 10 ))                                 # shear = 0.1 x compression, MPa
AMP=$(awk -v p=$P 'BEGIN{printf "%g", p/1000}') # GPa
SH=$(awk -v s=$S 'BEGIN{printf "%g", s/1000}')
TAG=lad_P${P}MPa_S${S}MPa
echo "task ${SLURM_ARRAY_TASK_ID:-0}: compression $AMP GPa, shear $SH GPa, tag $TAG"
bash ./run_comp.sh $AMP ${TEND:-10} 1 $TAG $SH
