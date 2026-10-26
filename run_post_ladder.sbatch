#!/bin/bash
#SBATCH -J ladpost
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --mem=12G
#SBATCH -t 04:00:00
#SBATCH --array=0-6
#SBATCH -o logs/ladpost_%A_%a.log
# Post-processing of the ladder levels (one array task per level, same order as run_ladder.sbatch):
# 14 animations (GIF + PNG frames), contact sheets, PNG zips and the data zip -> deliver_<tag>/
# Run after the ladder (or while it runs; it processes the frames that exist). Needs the python env from setup_py.sh.
set -uo pipefail
cd "$SLURM_SUBMIT_DIR"
module purge
module load python3 2>/dev/null || true
[ -f "$HOME/dynpost_venv/bin/activate" ] && source "$HOME/dynpost_venv/bin/activate"
python3 -c "import numpy, scipy, matplotlib, meshio, PIL" || { echo "FAILED: run 'bash setup_py.sh' once on the login node"; exit 1; }
LEVELS=(100 200 500 1000 2000 5000 10000)
P=${LEVELS[${SLURM_ARRAY_TASK_ID:-0}]}; S=$(( P / 10 )); TAG=lad_P${P}MPa_S${S}MPa
[ -s out_$TAG/$TAG.pvd ] && ls out_$TAG/${TAG}_s*.vtu >/dev/null 2>&1 || { echo "skip $TAG: no frames in out_$TAG"; exit 0; }
cd post
BASE=$TAG NP=2 python3 animate.py ../out_$TAG ../anim_$TAG
bash ./package.sh ../anim_$TAG ../deliver_$TAG ../out_$TAG ${TAG}_data.zip
