#!/bin/bash
# run.sh -- one restart-safe run.   usage: bash run.sh <speed m/s> [tend=30] [ramp=0.5] [hold=2.0] [vtkevery=10]
# env: MODE (mpi | serial, default mpi), tag.sh variables (OLD=1 or INERTIA=1: the set-up with rotational inertia j = 0.01), TSTOP (1400 K), RELAXCK, RELAXCHECK=0 (accept a relaxation checkpoint made with
#      other parameters), RELAXONLY=1 (used by run_relax.sh); MPI: NP, MPIRUN, PETSC / PETSCFALL / PETSCP1.
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
MODE=${MODE:-mpi}; case "$MODE" in mpi|serial) ;; *) echo "MODE must be mpi or serial"; exit 1;; esac
UPMS=${1:?plate speed in m/s}; TEND=${2:-30}; TR=${3:-0.5}; TH=${4:-2.0}; VE=${5:-10}
source ./tag.sh || exit 1
SC=$(awk -v d=$DT 'BEGIN{printf "%d", 0.01/d + 0.5}'); [ "$SC" -ge 1 ] || SC=1; VE=$((VE*SC)); CKE=$((50*SC))   # frames / checkpoints keep their spacing in ns
D=out_$TAG; B=$TAG
RUNKEY="$PHYS -dt $DT -pistonup $UN -pistonut $UT -ptr $TR -pth $TH"                  # what a resumed run must share with the run directory
ARGS="-pistonup $UN -pistonut $UT -ptr $TR -pth $TH -tend $TEND -dt $DT -vtkevery $VE -ckevery $CKE -tstop ${TSTOP:-1400} $PHYS"
# both set-ups run the v2 solver; each has its own relaxation checkpoint per RUNTAG (the relaxation depends on the parameters)
XS=_v2; RKS=$([ "$V2" = 1 ] && echo _v2 || echo _j)${RUNTAG:+_$RUNTAG}
# FASTASM=1 (default): tangent matrix assembled directly by the plugin on a fixed pattern (same matrix, results to round-off; FASTASM=0:
# FreeFEM forms). Exact slip law, local acceptance test (LOCALTOL, default 1e-10; 0 = off) and the lattice-rotation branch guard
# (BETAMAX radians, default pi/4; 0 = off). None of these is part of the settings key. VERBOSE=1: two timing lines per step.
FASTASM=${FASTASM:-1}; XARGS=""
if [ "$FASTASM" = 1 ]; then XARGS="-asm 1"; else XARGS="-asm 0 -fastasm 0"; fi   # FASTASM=0: mechanics tangent and staggered fields by the FreeFEM forms
XARGS="$XARGS -slipexact 1 -localtol ${LOCALTOL:-1e-10} -betamax ${BETAMAX:-0.785398}"
[ "${VERBOSE:-0}" = 1 ] && XARGS="$XARGS -verbose 1"
# ASMCHECK=N: every N-th Newton iteration of the run, the direct assembly is compared entry by entry with the form assembly (needs FASTASM=1)
if [ "${ASMCHECK:-0}" != 0 ]; then
  [ "$FASTASM" = 1 ] || { echo "FAILED: ASMCHECK needs FASTASM=1 (it checks the direct assembly against the forms)"; exit 1; }
  XARGS="$XARGS -asmcheck $ASMCHECK"
fi
RCK=${RELAXCK:-relax_checkpoint$RKS.txt}
# settings files written before the AD cleanup (2026-11) still carry -fdk/-fdq (finite-difference steps, no longer used): ignored
oldkey() { sed -E 's/ -fd[kq] [^ ]+//g' "$1"; }
keydiff() { diff <(oldkey "$1" 2>/dev/null | sed 's/ -/\n-/g') <(echo "$2" | sed 's/ -/\n-/g') | grep '^[<>]' | sed 's/^/    /'; }
rck_ok() { [ "${RELAXCHECK:-1}" = 0 ] || { [ -s $RCK.args ] && [ "$(oldkey $RCK.args)" = "$PHYS" ]; }; }
if [ "${RELAXONLY:-0}" = 1 ]; then
  # relaxation only: nothing to do if $RCK exists and was made with these parameters; otherwise run the -nrelax zero-load steps
  # from scratch at dt = 0.01 ns and publish the relaxed checkpoint (saved before the end-of-relaxation reset) as $RCK
  exec 9> $RCK.lock; flock 9   # one relaxation at a time; a second caller waits and then finds the file
  if [ -s $RCK ]; then
    rck_ok && { echo "relax checkpoint $RCK present (same parameters): nothing to do"; exit 0; }
    echo "FAILED: $RCK exists but was made with other parameters (< file, > now):"; keydiff $RCK.args "$PHYS"
    echo "  delete $RCK and $RCK.args to regenerate it, or set RELAXCHECK=0 to use it anyway"; exit 1
  fi
  TEND=0; D=out_relax$RKS; B=relax; rm -rf $D
  ARGS="-pistonup 0 -pistonut 0 -ptr $TR -pth $TH -tend 0 -dt 0.01 -vtkevery 0 -ckevery 1000000 -tstop ${TSTOP:-1400} $PHYS"
  echo "relaxation stage (no $RCK yet): $D/$B"
fi
mkdir -p $D
export FF_LOADPATH=${FF_LOADPATH:-$(dirname $(find $(dirname $(dirname $(which FreeFem++)))/lib* -name msh3.so 2>/dev/null | head -1) 2>/dev/null)}
done_p() { [ -s $D/${B}_log.csv ] && tail -1 $D/${B}_log.csv | cut -d, -f2 | awk -v e=$TEND -v d=$DT '{exit !($1+0==$1 && $1>=e-1.1*d)}'; }
done_p && { echo "already complete"; exit 0; }
[ -e $D/${B}_TSTOP.txt ] && { echo "ended at the temperature limit:"; cat $D/${B}_TSTOP.txt; exit 0; }
[ -e $D/${B}_NOCONV.txt ] && { echo "stopped earlier on a non-converged step:"; cat $D/${B}_NOCONV.txt; echo "(delete $D/${B}_NOCONV.txt to retry from its checkpoint, e.g. with a larger MAXITRETRY)"; exit 1; }
[ -e $D/${B}_NAN.txt ] && { echo "stopped earlier on NaN:"; cat $D/${B}_NAN.txt; echo "(last stage: $(sed -n 2p $D/${B}_stage.txt 2>/dev/null)); delete $D/${B}_NAN.txt to retry"; exit 1; }
[ -e $D/STALLED ] && { echo "STALLED earlier (see $D/STALLED); fix the cause, then delete $D/STALLED to retry"; exit 1; }
maxstep() { [ -s $D/${B}_log.csv ] && awk -F, '$1+0==$1 && $1>m {m=$1} END{print m+0}' $D/${B}_log.csv || echo 0; }
S0=$(maxstep)
# checkpoints are published through <ck>.manifest (immutable generations; the plain <ck> is only a convenience link)
CKF=$D/${B}_ck.txt
ckexists() { [ -s $CKF.manifest ] || [ -s $CKF ]; }
ckstep() {   # the committed checkpoint step, or fail (message on stderr) if the manifest/generation is unusable
  local g n m
  if [ -s $CKF.manifest ]; then
    n=$(sed -n 1p $CKF.manifest); g=$(sed -n 2p $CKF.manifest)
    [ -n "$g" ] && [ -s "$g" ] || { echo "ERROR: $CKF.manifest selects a missing/empty generation '$g'" >&2; return 1; }
    m=$(head -1 "$g")
  elif [ -s $CKF ]; then m=$(head -1 $CKF); n=$m
  else echo "ERROR: no checkpoint $CKF" >&2; return 1; fi
  case "$m" in ''|*[!0-9]*) echo "ERROR: checkpoint step '$m' is not an integer" >&2; return 1;; esac
  [ "$n" = "$m" ] || { echo "ERROR: manifest step '$n' differs from the generation's step '$m'" >&2; return 1; }
  echo "$m"
}
if [ "${RELAXONLY:-0}" != 1 ]; then
  # settings of this run directory: written at the first start, a resume with other settings is refused
  if [ -s $D/$B.args ]; then
    [ "$(oldkey $D/$B.args)" = "$RUNKEY" ] || { echo "FAILED: $D was started with other settings (< there, > now):"; keydiff $D/$B.args "$RUNKEY"; echo "  use another RUNTAG for a run with other settings"; exit 1; }
  elif ckexists; then echo "FAILED: $D has a checkpoint but no $B.args (not made by these scripts): use another RUNTAG"; exit 1
  else
    [ -s $RCK ] || { echo "FAILED: $RCK missing: run 'MODE=$MODE bash run_relax.sh' first (the submit scripts schedule it automatically)"; exit 1; }
    rck_ok || { echo "FAILED: $RCK was made with other parameters than this run (< file, > now):"; keydiff $RCK.args "$PHYS"; echo "  regenerate it (delete $RCK $RCK.args, then bash run_relax.sh) or set RELAXCHECK=0"; exit 1; }
    cp $RCK $D/${B}_ck.txt; echo "$RUNKEY" > $D/$B.args
  fi
  # solver options outside the settings key (assembly, slip-law evaluation, acceptance tests, output): one line per start in
  # <tag>.xargs (UTC time | options), so every continuation is traceable; a change from the previous start is reported
  if [ -s $D/$B.xargs ] && [ "$(tail -n 1 $D/$B.xargs | cut -d'|' -f2-)" != "$XARGS" ]; then
    echo "note: solver options outside the settings key differ from the previous start: was '$(tail -n 1 $D/$B.xargs | cut -d'|' -f2-)', now '$XARGS' ($D/$B.xargs)"
  fi
  echo "$(date -u +%FT%TZ)|$XARGS" >> $D/$B.xargs
fi
R=""; ckexists && R="-restart $CKF"
if [ -n "$R" ]; then   # resuming: drop log/diag/thermo/power/newton rows written after the checkpoint (they are recomputed, or belonged to an abandoned branch), so every file holds the active history once
  NCK=$(ckstep) || { echo "FAILED: unusable checkpoint; logs left untouched (run recover_ck.sh or fix the manifest)"; exit 1; }
  for f in $D/${B}_log.csv $D/${B}_diag.csv $D/${B}_thermo.csv $D/${B}_power.csv $D/${B}_newton.csv; do
    [ -s "$f" ] && awk -F, -v n=$NCK 'NR==1 && $1=="step" {print; next} $1+0<=n' "$f" > "$f.tmp" && mv "$f.tmp" "$f"
  done
  echo "resuming from the checkpoint at step $NCK (furthest step reached so far: $S0)"
  HDR="step,t_ns,load_GPa,iters,locres,KE_trans,KE_rot,E_elastic,E_stored,Q_heat,W_ext,D_cum,audit_residual,max_abs_slip,mean_abs_slip,max_dT,max_abs_dtheta_deg,mean_rhoS_um2,max_rhoG_um2,min_phi,neg_D_count,W_plate,audit_full"
  [ -s $D/${B}_log.csv ] || echo "$HDR" > $D/${B}_log.csv   # (the solver writes the header only on a fresh start)
fi
touch $D/$B.pvd
# the plugin kwc2d.so (element-map derivatives by AD) and cosserat_kwc_mpi.edp, rebuilt only when stale
bash plugin/build_plugin.sh >> $D/$B.out 2>&1 || { echo "FAILED: plugin/build_plugin.sh (see plugin/build.log and $D/$B.out)"; exit 1; }
# the hardware this link runs on (Longleaf nodes differ; 'grep HOST out_*/*.out' compares runs): node, CPU model, the CPUs Slurm gave
# this job and how many physical cores they are (fewer cores than ranks = ranks share cores through hyperthreading), node load
hostinfo() {
  local cpus model phys
  cpus=$(taskset -cp $$ 2>/dev/null | sed 's/.*: //')
  model=$(lscpu 2>/dev/null | sed -n 's/^Model name: *//p' | head -1)
  phys=$(lscpu -p=CPU,CORE,SOCKET 2>/dev/null | grep -v '^#' | awk -F, -v L="$cpus" '
    BEGIN { n = split(L, a, ","); for (i = 1; i <= n; i++) { if (split(a[i], r, "-") == 2) { for (c = r[1]; c <= r[2]; c++) ok[c] = 1 } else ok[a[i]] = 1 } }
    ($1 in ok) && !(($3 ":" $2) in u) { u[$3 ":" $2] = 1; k++ } END { print k + 0 }')
  HPHYS=${phys:-0}
  echo "[$(date +%T)] HOST $(hostname -s) | ${model:-?} | job ${SLURM_JOB_ID:-none} cpus ${cpus:-?} = ${phys:-?} physical cores for ${NP:-1} ranks | load $(cut -d' ' -f1-3 /proc/loadavg) | step $S0"
}
[ "$MODE" = mpi ] && NP=${NP:-${SLURM_NTASKS:-4}}
hostinfo >> $D/$B.out
echo "[$(date +%T)] code version $(cat VERSION 2>/dev/null || echo unknown)" >> $D/$B.out
# frozen solver (version 8.0): the solver sources, the generated files and tag.sh against SOLVER_SHA256 (checks/verify_solver.sh)
if [ -s SOLVER_SHA256 ]; then
  if sha256sum -c --quiet SOLVER_SHA256 >/dev/null 2>&1; then echo "[$(date +%T)] solver frozen: files match SOLVER_SHA256" >> $D/$B.out
  else echo "[$(date +%T)] solver MODIFIED: $(sha256sum -c --quiet SOLVER_SHA256 2>/dev/null | tr '\n' ' ')" >> $D/$B.out; fi
fi
# hardware guard: every rank needs its own physical core. --hint=nomultithread does not always hold, and ranks sharing cores through
# hyperthreading ran about 10x slower. Stops at once with the resubmit line; HWGUARD=0 runs anyway. (Restart-safe: nothing is lost.)
if [ "$MODE" = mpi ] && [ "${HWGUARD:-1}" = 1 ] && [ "${HPHYS:-0}" -gt 0 ] && [ "$HPHYS" -lt "${NP:-1}" ]; then
  echo "[$(date +%T)] HARDWARE: $HPHYS physical cores for $NP ranks on $(hostname -s) (hyperthreading); not running. Resubmit with HWEXCLUDE=$(hostname -s) (and HWCONSTRAINT / HWEXCL=1, see hw.sh), or HWGUARD=0 to run anyway" | tee -a $D/$B.out
  exit 3
fi
if [ "$MODE" = mpi ]; then
  EDP=cosserat_kwc${XS}_mpi.edp
  NP=${NP:-${SLURM_NTASKS:-4}}
  if [ -z "${MPIRUN:-}" ]; then
    # prefer the mpiexec FreeFem++-mpi was built with (recorded in ff-mpirun); a foreign launcher (e.g. srun without PMI)
    # starts NP independent single-rank copies instead of one NP-rank job (the solver stops on -npexpect)
    FFM=$(command -v ff-mpirun 2>/dev/null); ME=""
    [ -n "$FFM" ] && ME=$(grep -o 'mpirun="[^"]*"' "$FFM" | tail -1 | cut -d'"' -f2)
    # --bind-to none: Slurm may hand out hyperthreads of fewer physical cores than ranks, and OpenMPI refuses to bind then
    if [ -n "$ME" ] && [ -x "$ME" ]; then MPIRUN="$ME -n $NP --bind-to none"
    elif command -v mpiexec >/dev/null; then MPIRUN="mpiexec -n $NP --bind-to none"
    elif command -v mpirun >/dev/null; then MPIRUN="mpirun -np $NP"
    elif [ -n "${SLURM_JOB_ID:-}" ] && command -v srun >/dev/null; then MPIRUN="srun -n $NP"
    else echo "FAILED: no MPI launcher found; load the MPI module FreeFem++-mpi was built with, or set MPIRUN"; exit 1; fi
  fi
  command -v FreeFem++-mpi >/dev/null || { echo "FAILED: FreeFem++-mpi not on PATH"; exit 1; }
  PARGS=(); [ -n "${PETSC:-}" ] && PARGS+=(-petsc "$PETSC"); [ -n "${PETSCFALL:-}" ] && PARGS+=(-petscfall "$PETSCFALL"); [ -n "${PETSCP1:-}" ] && PARGS+=(-petscp1 "$PETSCP1"); [ -n "${PETSCPH:-}" ] && PARGS+=(-petscph "$PETSCPH")
  echo "[$(date +%T)] $MPIRUN FreeFem++-mpi $EDP ($NP ranks)" >> $D/$B.out
  L0=$(wc -l < $D/$B.out)
  $MPIRUN FreeFem++-mpi -nw -v 0 $EDP $ARGS $XARGS "${PARGS[@]}" -npexpect $NP -out $D/$B $R >> $D/$B.out 2>&1
  EXIT=$?
  tail -n +$((L0+1)) $D/$B.out | grep -q "LAUNCHER ERROR" && { echo "launcher mismatch (see $D/$B.out): not a solver failure; set MPIRUN and resubmit"; grep -m1 "LAUNCHER ERROR" $D/$B.out; exit 1; }
else
  EDP=cosserat_kwc$XS.edp
  if ! ls ${FF_LOADPATH:-/nonexistent}/IncompleteCholesky.so >/dev/null 2>&1 && ! find $(dirname $(dirname $(which FreeFem++)))/lib* -name IncompleteCholesky.so 2>/dev/null | grep -q .; then
    echo "FAILED: the FreeFem plugin IncompleteCholesky (ILU preconditioner of the serial GMRES solve) was not found"; exit 1; fi
  echo "[$(date +%T)] FreeFem++ $EDP (serial)" >> $D/$B.out
  FreeFem++ -nw -v 0 $EDP $ARGS $XARGS -out $D/$B $R >> $D/$B.out 2>&1
  EXIT=$?
fi
[ "$EXIT" != 0 ] && [ -s $D/${B}_stage.txt ] && { echo "solver exit code $EXIT; last stage reached:"; cat $D/${B}_stage.txt; }
S1=$(maxstep)
if done_p; then echo COMPLETE
elif [ -e $D/${B}_TSTOP.txt ]; then echo "ENDED at the temperature limit:"; cat $D/${B}_TSTOP.txt
elif [ -e $D/${B}_NAN.txt ]; then echo "STOPPED on NaN:"; cat $D/${B}_NAN.txt
elif [ -e $D/${B}_NOCONV.txt ]; then echo "STOPPED on a non-converged step:"; cat $D/${B}_NOCONV.txt
elif [ "$S1" -le "$S0" ]; then
  { echo "STALLED: no new steps beyond step $S0 (solver exit code $EXIT) at $(date)"; echo "--- last lines of $D/$B.out ---"; tail -n 25 $D/$B.out; } | tee $D/STALLED
  exit 1
else echo "INCOMPLETE at step $S1 (exit code $EXIT): resubmit (restart-safe)"; fi
if [ "${RELAXONLY:-0}" = 1 ]; then
  [ -s $D/${B}_relaxed_ck.txt ] || { echo "FAILED: the relaxation did not write $D/${B}_relaxed_ck.txt (see $D/$B.out)"; exit 1; }
  echo "$PHYS" > $RCK.args.tmp && cp $D/${B}_relaxed_ck.txt $RCK.tmp && mv $RCK.args.tmp $RCK.args && mv $RCK.tmp $RCK && echo "relax checkpoint written: $RCK (parameters in $RCK.args)"
  echo "$(date -u +%FT%TZ)|$XARGS" > $RCK.xargs   # options outside the key used by the relaxation (record only)
fi
