#!/bin/bash
#SBATCH -J grainSum
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=4G
#SBATCH -t 00:30:00
#SBATCH -o logs/grainSum_%j.log
# Cross-run grain summary (refinement factor, grain count, mean size +- std, VL/L/H misorientation) over angles x speeds, from the
# deliver_<tag>/grains/ folders of the post arrays (missing runs are skipped) -> grain_summary<suffix>/.
# env: ANGLES ("0 45 90"), UPS ("25 50 100 150 250 400"), RUNTAG, DT, MODE (mpi | serial) as for the runs; MISOR (sym|raw), AREA (territory|plateau), SMOOTH, TAVG
set -uo pipefail
cd "$SLURM_SUBMIT_DIR"
module purge; module load python3 2>/dev/null || true
[ -f "$HOME/dynpost_venv/bin/activate" ] && source "$HOME/dynpost_venv/bin/activate"
python3 -c "import numpy, pandas, matplotlib" || { echo "FAILED: source ~/dynpost_venv/bin/activate; pip install pandas"; exit 1; }
MODE=${MODE:-mpi}; UPMS=""; source ./tag.sh || exit 1
python3 post/grain_summary.py --root . --prefix deliver --suffix "$SUFFIX" --angles ${ANGLES:-0 45 90} --ups ${UPS:-25 50 100 150 250 400} \
  --misor ${MISOR:-sym} --area ${AREA:-territory} --smooth ${SMOOTH:-1} --tavg ${TAVG:-0} --out grain_summary$SUFFIX
