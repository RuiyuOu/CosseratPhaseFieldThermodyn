#!/bin/bash
#SBATCH -J postMPI
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --mem=12G
#SBATCH -t 06:00:00
#SBATCH --array=0-5
#SBATCH -o logs/postMPI_%A_%a.log
# Post-processing of the MPI runs (tag suffix _mpi; SHARETAG=1: serial tag) of the mobile-drag piston ladder (tasks 0-5, variant A). Set the same DRAGRHO/VC/GND/RECOV/ETAPHIMIN/ETASTAR (and CODE/MUCDUAL/DUALR/KWCEDIFF/KWCFIELD/KWCONEWAY) as for the run (default 3/2/curl/2/1e-2/1 -> out_pistM3v2GRcV_gc_ne_pm1e-2_es1_*) -> deliver_<tag>/
set -uo pipefail
cd "$SLURM_SUBMIT_DIR"
module purge
module load python3 2>/dev/null || true
[ -f "$HOME/dynpost_venv/bin/activate" ] && source "$HOME/dynpost_venv/bin/activate"
python3 -c "import numpy, scipy, matplotlib, meshio, PIL" || { echo "FAILED: run 'bash setup_py.sh' once on the login node"; exit 1; }
UPS=(25 50 100 150 250 400)
if [ -n "${GNDCURL:-}" ] || [ -n "${GNDPHI2:-}" ]; then GP=${GNDPHI2:-1}; GCL=${GNDCURL:-0}   # legacy flags win
else GND=${GND:-curl}; case $GND in curl) GP=0; GCL=1;; grad) GP=0; GCL=0;; phigrad) GP=1; GCL=0;; *) echo "GND must be curl|grad|phigrad" >&2; exit 1;; esac; fi
ID=${SLURM_ARRAY_TASK_ID:-0}; SFX="M${DRAGRHO:-3}v${VC:-2}"; [ "$GP" = 0 ] && SFX=${SFX}G; [ "${RECOV:-2}" != 0 ] && SFX=${SFX}R; { [ "${DISSIP:-1}" != 0 ] || [ "${CDSTORE:-1}" != 0 ]; } && SFX=${SFX}c; case ${EOS:-vinet} in vinet) SFX=${SFX}V;; murnaghan) SFX=${SFX}Mu;; esac; { [ "${DT:-0.01}" != 0.01 ] || [ "${NREF:-20}" != 20 ]; } && SFX=${SFX}_dt${DT:-0.01}_nr${NREF:-20}; [ "${SPINMAX:-0}" != 0 ] && SFX=${SFX}_sm${SPINMAX}; [ "${SYMFOLD:-0}" != 0 ] && SFX=${SFX}_sf${SYMFOLD}; [ "$GCL" != 0 ] && SFX=${SFX}_gc; [ "${CODE:-stable}" = exp ] && SFX=${SFX}_exp; [ "${MUCDUAL:-0}" != 0 ] && SFX=${SFX}_du; [ "${MUCDUAL:-0}" != 0 ] && [ "${DUALR:-0}" != 0 ] && SFX=${SFX}_dr${DUALR}; [ "${KWCEDIFF:-0}" = 0 ] && SFX=${SFX}_ne; [ "${KWCE:-0.004}" != 0.004 ] && SFX=${SFX}_e${KWCE}; [ "${KWCREG:-0}" != 0 ] && SFX=${SFX}_kg; [ "${KWCREG:-0}" != 0 ] && [ "${KWCGAM:-0}" != 0 ] && SFX=${SFX}${KWCGAM}; [ "${KWCFIELD:-0}" != 0 ] && SFX=${SFX}_kf; [ "${KWCONEWAY:-0}" != 0 ] && SFX=${SFX}_ow; PMN=${ETAPHIMIN:-1e-2}; [ "$PMN" != 0 ] && SFX=${SFX}_pm${PMN}; [ "${QPHI:-100e3}" != 100e3 ] && SFX=${SFX}_q${QPHI}; [ "${ETAPHI:-0.1}" != 0.1 ] && SFX=${SFX}_ep${ETAPHI}; EST=${ETASTAR:-1}; [ "$EST" != 1e-3 ] && SFX=${SFX}_es${EST}; [ -n "${RUNTAG:-}" ] && SFX=${SFX}_${RUNTAG}; [ "${SHARETAG:-0}" = 1 ] || SFX=${SFX}_mpi; PFX=pist$SFX
ATAG=""; [ "${ANGLE:-45}" != confined ] && ATAG=_a${ANGLE:-45}
UP=${UPS[$ID]}; TAG=${PFX}${ATAG}_up${UP}ms
[ -s out_$TAG/$TAG.pvd ] && ls out_$TAG/${TAG}_s*.vtu >/dev/null 2>&1 || { echo "skip $TAG: no frames in out_$TAG"; exit 0; }
cd post
BASE=$TAG NP=2 python3 animate.py ../out_$TAG ../anim_$TAG
# animations: incremental (frames already rendered in ../anim_$TAG are kept; FORCE=1 redoes all, RELIM=1 also recomputes the colour limits)
# grain count + boundary misorientation history (../grains_$TAG) and the grain_theta / grain_territory / grain_misorientation animations (into ../anim_$TAG,
# packaged like the others): incremental too (GRAD, MERGE, MINAREA, SYM, LO, MID, HI env)
BASE=$TAG NP=2 python3 grain_post.py ../out_$TAG ../grains_$TAG ../anim_$TAG || echo "grain_post failed (animations are unaffected)"
bash ./package.sh ../anim_$TAG ../deliver_$TAG ../out_$TAG ${TAG}_data.zip ../grains_$TAG
