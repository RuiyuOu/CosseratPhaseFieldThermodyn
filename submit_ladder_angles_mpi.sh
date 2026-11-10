#!/bin/bash
# Full variant-A velocity ladder (25 50 100 150 250 400 m/s) at plate angles 0, 45 and 90 degrees, MPI/PETSc, 16 ranks and 48 GB per run, with the current
# defaults (3-point Gauss rule on the finite-strain energy -> tag "Q3", exact tangent, KKT/PDAS phase field, -safe 1, DT=0.01, TEND=30), i.e. the same settings as
#     pistM3v2GRcVQ3_gc_ne_pm1e-2_es1_mpi_a0_up400ms
# and the post-processing (animations, grain count, tau/T/pressure/mismatch statistics, packaging) scheduled behind each angle's last run link.
# Outputs: out_pistM3v2GRcVQ3_gc_ne_pm1e-2_es1_mpi_a<ANGLE>_up<speed>ms/   post: anim_<tag>/ grains_<tag>/ tau_<tag>/ deliver_<tag>/
#
# usage (from this directory on Longleaf):   bash submit_ladder_angles_mpi.sh [N_chained_24h_arrays=10]
#   env: ANGLES ("0 45 90"), ARRAY ("0-5" = 25 50 100 150 250 400 m/s; e.g. ARRAY=3-5 for 150 250 400), NP (16), MEM (48G), TEND (30), SAFE (1),
#        RUNTAG (empty; set to keep another series apart), SERIAL_ANGLES=1 (start the next angle only after the previous angle's post-processing has finished,
#        so at most 6 x 16 cores are in use at a time; default 0 = all three angles are queued together, 3 x 6 x 16 = 288 cores)
#        post-processing: TAUWIN0 TAUWIN1 TAUMARK TAUEVERY (see run_post_piston_mob_mpi.sbatch)
# grain summary figures (post/grain_summary.py) are scheduled after all post arrays: run_grain_summary.sbatch
# prerequisites: cosserat2d_mob_stable_mpi.edp, out_A/A_dyn_relaxed_ck.txt, the MPI-enabled freefem module used by run_piston_mob_mpi.sbatch, 'bash setup_py.sh' run once (post).
# Every level is restart-safe: running this script again resumes unfinished levels and skips finished ones (a complete level exits at once).
set -euo pipefail
cd "$(dirname "$0")"
N=${1:-10}
export NP=${NP:-16} MEM=${MEM:-48G} ARRAY=${ARRAY:-0-5} SAFE=${SAFE:-1}
[ -s cosserat2d_mob_stable_mpi.edp ] || python3 mpi/make_mpi.py cosserat2d_mob_stable.edp cosserat2d_mob_stable_mpi.edp
[ -s out_A/A_dyn_relaxed_ck.txt ] || { echo "missing out_A/A_dyn_relaxed_ck.txt"; exit 1; }
mkdir -p logs
dep=""; alldeps=""
for A in ${ANGLES:-0 45 90}; do
  echo "== angle $A: ARRAY=$ARRAY, $NP ranks, $MEM"
  if [ "${SERIAL_ANGLES:-0}" = 1 ] && [ -n "$dep" ]; then export SBATCH_DEPENDENCY="afterany:$dep"; fi
  out=$(ANGLE=$A bash run_all_A_mpi.sh "$N")
  echo "$out"
  dep=$(echo "$out" | sed -n 's/^post-processing array: \([0-9]*\).*/\1/p')
  [ -n "$dep" ] && alldeps="$alldeps:$dep"
done
[ -n "$alldeps" ] && echo "grain summary: $(sbatch --parsable --export=ALL --dependency=afterany$alldeps run_grain_summary.sbatch)  -> grain_summary_*/ (refinement, count, size, VL/L/H misorientation)"
echo "progress:  tail -qn 1 out_pistM3v2GRcVQ3_gc_ne_pm1e-2_es1_mpi_a*_up*ms/*_log.csv     queue: squeue -u \$USER"
