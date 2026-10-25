#!/bin/bash
#SBATCH -J dynF
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=8G
#SBATCH -t 12:00:00
#SBATCH -o logs/dynF_%j.log
set -euo pipefail
RUN=run_dyn.sh
for d in "$SLURM_SUBMIT_DIR" "$(dirname "$(scontrol show job "$SLURM_JOB_ID" 2>/dev/null | sed -n 's/^ *Command=//p' | head -1)")"; do
  [ -f "$d/$RUN" ] && { cd "$d"; break; }
done
[ -f "$RUN" ] || { echo "FAILED: cannot find $RUN; submit from inside the bundle folder"; exit 1; }
module purge
module use "$HOME/modulefiles"
module load freefem/4.15 sqlite/3.45.0
which FreeFem++ g++
mkdir -p logs
bash ./run_dyn.sh F
