#!/bin/bash
#SBATCH -J postM
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --mem=12G
#SBATCH -t 06:00:00
#SBATCH --array=0-5
#SBATCH -o logs/postM_%A_%a.log
# Post-processing of the mobile-drag piston ladder (tasks 0-5, variant A). Set the same DRAGRHO/VC/GNDPHI2/RECOV as for the run (default 3/2/1/2 -> out_pistM3v2R_up*) -> deliver_<tag>/
set -uo pipefail
cd "$SLURM_SUBMIT_DIR"
module purge
module load python3 2>/dev/null || true
[ -f "$HOME/dynpost_venv/bin/activate" ] && source "$HOME/dynpost_venv/bin/activate"
python3 -c "import numpy, scipy, matplotlib, meshio, PIL" || { echo "FAILED: run 'bash setup_py.sh' once on the login node"; exit 1; }
UPS=(25 50 100 150 250 400)
ID=${SLURM_ARRAY_TASK_ID:-0}; SFX="M${DRAGRHO:-3}v${VC:-2}"; [ "${GNDPHI2:-1}" = 0 ] && SFX=${SFX}G; [ "${RECOV:-2}" != 0 ] && SFX=${SFX}R; { [ "${DISSIP:-1}" != 0 ] || [ "${CDSTORE:-1}" != 0 ]; } && SFX=${SFX}c; case ${EOS:-vinet} in vinet) SFX=${SFX}V;; murnaghan) SFX=${SFX}Mu;; esac; { [ "${DT:-0.01}" != 0.01 ] || [ "${NREF:-20}" != 20 ]; } && SFX=${SFX}_dt${DT:-0.01}_nr${NREF:-20}; [ "${SPINMAX:-0}" != 0 ] && SFX=${SFX}_sm${SPINMAX}; PFX=pist$SFX
UP=${UPS[$ID]}; TAG=${PFX}_up${UP}ms
[ -s out_$TAG/$TAG.pvd ] && ls out_$TAG/${TAG}_s*.vtu >/dev/null 2>&1 || { echo "skip $TAG: no frames in out_$TAG"; exit 0; }
cd post
BASE=$TAG NP=2 python3 animate.py ../out_$TAG ../anim_$TAG
bash ./package.sh ../anim_$TAG ../deliver_$TAG ../out_$TAG ${TAG}_data.zip
