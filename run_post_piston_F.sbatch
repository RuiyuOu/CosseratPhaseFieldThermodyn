#!/bin/bash
#SBATCH -J pistFpost
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --mem=12G
#SBATCH -t 06:00:00
#SBATCH --array=0-5
#SBATCH -o logs/pistFpost_%A_%a.log
# Post-processing of the piston ladder levels (variant F only, same order as run_piston_F.sbatch): 14 animations, contact sheets, PNG zips, data zip -> deliver_<tag>/
set -uo pipefail
cd "$SLURM_SUBMIT_DIR"
module purge
module load python3 2>/dev/null || true
[ -f "$HOME/dynpost_venv/bin/activate" ] && source "$HOME/dynpost_venv/bin/activate"
python3 -c "import numpy, scipy, matplotlib, meshio, PIL" || { echo "FAILED: run 'bash setup_py.sh' once on the login node"; exit 1; }
UPS=(25 50 100 150 250 400)
ID=${SLURM_ARRAY_TASK_ID:-0}; PFX=pistF
UP=${UPS[$ID]}; TAG=${PFX}_up${UP}ms
[ -s out_$TAG/$TAG.pvd ] && ls out_$TAG/${TAG}_s*.vtu >/dev/null 2>&1 || { echo "skip $TAG: no frames in out_$TAG"; exit 0; }
cd post
BASE=$TAG NP=2 python3 animate.py ../out_$TAG ../anim_$TAG
bash ./package.sh ../anim_$TAG ../deliver_$TAG ../out_$TAG ${TAG}_data.zip
