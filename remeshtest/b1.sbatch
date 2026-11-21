#!/bin/bash
#SBATCH -J b1curv
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=8G
#SBATCH -t 04:00:00
#SBATCH -o logs/b1_%j.log
# B1: lattice-curvature refinement indicator on the VTU frames of 2D runs (remesh/b1_curvature.py), one run after the other.
# usage (repo root):  sbatch remesh/b1.sbatch out_a90_up400ms out_a0_up400ms ...      (works on runs still in progress: frames so far)
# env: HMIN (0.1), HMAX (0.4), GRADE (0.3), EVERY (1 = every frame), B1ARGS (more options, e.g. "--dstar 0.5 1 --buffer 0.5")
# output: remesh/b1_<tag>/b1_frames.csv, b1_summary.png, b1_hfinal_*.png
set -uo pipefail
cd "$SLURM_SUBMIT_DIR"
module purge; module load python3 2>/dev/null || true
[ -f "$HOME/dynpost_venv/bin/activate" ] && source "$HOME/dynpost_venv/bin/activate"
python3 -c "import numpy, scipy, matplotlib, meshio" || { echo "FAILED: run 'bash setup_py.sh' once on the login node"; exit 1; }
[ $# -ge 1 ] || { echo "give one or more run folders (out_<tag>)"; exit 1; }
for R in "$@"; do
  T=$(basename "$R"); T=${T#out_}
  echo "== $R"
  python3 remesh/b1_curvature.py "$R" --hmin ${HMIN:-0.1} --hmax ${HMAX:-0.4} --grade ${GRADE:-0.3} --every ${EVERY:-1} --out remesh/b1_$T ${B1ARGS:-}
done
