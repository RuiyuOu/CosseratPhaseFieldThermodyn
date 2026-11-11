#!/bin/bash
#SBATCH -J post
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --mem=12G
#SBATCH -t 06:00:00
#SBATCH --array=0-5
#SBATCH -o logs/post_%A_%a.log
# Post-processing, one speed per array task (0-5 = 25 50 100 150 250 400 m/s): animations (anim_<tag>/) and the delivery folder
# deliver_<tag>/ (contact sheets, gifs, png zips, data zip, grains/ = grain history, tau/ = stress statistics).  Set the same MODE, ANGLE, RUNTAG, DT (SHARETAG) as for the runs.
# env: TAUWIN0 TAUWIN1 (loading window, default 0.5-30 ns), TAUEVERY (every N-th snapshot);
#      grain_post: GRAD, MERGE, MINAREA, SYM, LO, MID, HI; animations: FORCE=1 redoes all frames, RELIM=1 also the colour limits.
set -uo pipefail
cd "$SLURM_SUBMIT_DIR"
module purge
module load python3 2>/dev/null || true
[ -f "$HOME/dynpost_venv/bin/activate" ] && source "$HOME/dynpost_venv/bin/activate"
python3 -c "import numpy, scipy, matplotlib, meshio, PIL" || { echo "FAILED: run 'bash setup_py.sh' once on the login node"; exit 1; }
UPS=(25 50 100 150 250 400)
ID=${SLURM_ARRAY_TASK_ID:-0}
MODE=${MODE:-mpi}; UPMS=${UPS[$ID]}; source ./tag.sh || exit 1
[ -s out_$TAG/$TAG.pvd ] && ls out_$TAG/${TAG}_s*.vtu >/dev/null 2>&1 || { echo "skip $TAG: no frames in out_$TAG"; exit 0; }
cd post
BASE=$TAG NP=2 python3 animate.py ../out_$TAG ../anim_$TAG
mkdir -p ../deliver_$TAG
BASE=$TAG NP=2 python3 grain_post.py ../out_$TAG ../deliver_$TAG/grains ../anim_$TAG || echo "grain_post failed (animations are unaffected)"
LOGARG=""; [ -s ../out_$TAG/${TAG}_log.csv ] && LOGARG="--log ../out_$TAG/${TAG}_log.csv"
python3 tau_stats.py ../out_$TAG --out ../deliver_$TAG/tau --dt $DT --window ${TAUWIN0:-0.5} ${TAUWIN1:-30} --every ${TAUEVERY:-1} $LOGARG || echo "tau_stats failed (animations are unaffected)"
bash ./package.sh ../anim_$TAG ../deliver_$TAG ../out_$TAG ${TAG}_data.zip
