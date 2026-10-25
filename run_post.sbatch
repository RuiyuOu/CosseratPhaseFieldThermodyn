#!/bin/bash
#SBATCH -J dynpost
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=16G
#SBATCH -t 04:00:00
#SBATCH -o logs/dynpost_%j.log
set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
module purge
module load python3 2>/dev/null || true
[ -f "$HOME/dynpost_venv/bin/activate" ] && source "$HOME/dynpost_venv/bin/activate"   # created by setup_py.sh
python3 -c "import numpy, scipy, matplotlib, meshio, PIL" || { echo "FAILED: run 'bash setup_py.sh' once on the login node"; exit 1; }
cd post
for V in A F; do
  [ -d ../out_$V ] || { echo "skip $V: no out_$V"; continue; }
  BASE=${V}_dyn NP=4 python3 animate.py ../out_$V ../anim_$V
  bash ./package.sh ../anim_$V ../deliver_$V ../out_$V ${V}_dyn_data.zip
done
