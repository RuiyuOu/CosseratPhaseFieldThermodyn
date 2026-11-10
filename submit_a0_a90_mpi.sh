#!/bin/bash
# Submit the patched solver (Gauss-point elastic energy, exact tangent, KKT/PDAS phase field, -safe 1) for the 400 m/s plate impact at ANGLE = 0 and 90 degrees,
# MPI/PETSc, 16 ranks per run, everything else exactly as the earlier production run  pistM3v2GRcV_gc_ne_pm1e-2_es1_mpi_a90_up400ms
# (DRAGRHO=3 VC=2 GND=curl RECOV=2 EOS=vinet DISSIP=1 CDSTORE=1 KWCEDIFF=0 ETAPHIMIN=1e-2 ETASTAR=1, DT=0.01, TEND=30; the launcher defaults).
# New runs get their own directories (the tag gains "Q3" for the 3-point Gauss rule), so the old output is never touched or resumed:
#     out_pistM3v2GRcVQ3_gc_ne_pm1e-2_es1_mpi_a0_up400ms      out_pistM3v2GRcVQ3_gc_ne_pm1e-2_es1_mpi_a90_up400ms
# usage (from this directory on Longleaf):   bash submit_a0_a90_mpi.sh [N_chained_24h_arrays=10]
#   env: NP (16), MEM (48G), TEND (30), ANGLES ("0 90"), RUNTAG (empty; set e.g. RUNTAG=k1 to keep another series apart), SAFE (1)
# prerequisites: cosserat2d_mob_stable_mpi.edp (generated from the serial file: python3 mpi/make_mpi.py cosserat2d_mob_stable.edp cosserat2d_mob_stable_mpi.edp)
#                out_A/A_dyn_relaxed_ck.txt (the relaxed start, as before), the MPI-enabled freefem module used by run_piston_mob_mpi.sbatch.
# Then compare with the old run:   bash compare_a0_a90.sh
set -euo pipefail
cd "$(dirname "$0")"
N=${1:-10}; export NP=${NP:-16} MEM=${MEM:-48G} ARRAY=5 SAFE=${SAFE:-1}
[ -s cosserat2d_mob_stable_mpi.edp ] || python3 mpi/make_mpi.py cosserat2d_mob_stable.edp cosserat2d_mob_stable_mpi.edp
[ -s out_A/A_dyn_relaxed_ck.txt ] || { echo "missing out_A/A_dyn_relaxed_ck.txt"; exit 1; }
for A in ${ANGLES:-0 90}; do
  echo "== angle $A: 400 m/s, $NP ranks"
  ANGLE=$A bash run_all_A_mpi.sh "$N"
done
echo "check progress:  tail -n 1 out_pistM3v2GRcVQ3_gc_ne_pm1e-2_es1_mpi_a{0,90}_up400ms/*_log.csv"
