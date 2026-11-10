#!/bin/bash
#   ANGLE (default 45): plate-impact angle in degrees from -x. The plate speed V (first argument) splits into u_n = V cos(ANGLE) (compression, -x) and
#     u_t = V sin(ANGLE) (shear, +y): 0 = pure compression, 90 = pure shear. Plate BCs (-dir shearpiston): left edge clamped, right edge the plate,
#     top/bottom traction-free. Tag gets _a<ANGLE>. ANGLE=confined = the old confined piston (-dir piston: rollers top/bottom, frictionless wall), old tag.
# ===== MPI/PETSc version: same options, tags and file formats as run_piston_mob.sh, solver cosserat2d_mob_mpi.edp =====
#   NP (default $SLURM_NTASKS, else 4) MPI ranks; MPIRUN launcher (default: the mpiexec recorded in ff-mpirun, else mpiexec/mpirun on PATH, else srun);
#   PETSC: KSP/PC options for the (u, theta) Newton system (default GMRES + additive Schwarz with exact subdomain LU);
#   PETSCFALL: used for a solve that misses LINTOL (default MUMPS direct); PETSCP1: phase-field/heat systems.
#   The tag gets the suffix _mpi so the serial runs are never touched; SHARETAG=1 drops it, so an MPI job resumes
#   (and continues) the serial run of the same tag -- the checkpoints are identical in format.
#   LINSOLVER/PRECEVERY are ignored (PETSc solves the systems).
# Mobile-density drag variant (cosserat2d_mob.edp). Env: DRAGRHO (default 3 = evolved mobile density), VC (default 2 km/s relativistic drag; 0 = off),
#   RHOMSAT (100 um^-2), C2M (10), RHOFCUT (1e4 um^-2), RHOMMIN (1), GNDPHI2 (default 1 = original), RECOV (default 2 = Ask 2019 Eq. 80: Heaviside(phi_dot) tanh(C_A^2|grad theta|^2); 0 = original |dphi|), CA (sqrt 10 um), CD (default 100 for RECOV=2, else 10), RELAXFIX (default 1), DT (default 0.01 ns), NREF (default 20 steps),
#   BCKIN (default 1: piston face/rollers/wall carry their exact prescribed velocity and acceleration instead of the Newmark recursion, which grows linearly there; 0 = old),
#   PREDICTOR (default vel: Newton starts from U + dt V; avg = old U + dt V + dt^2/4 A; none), SPINBRANCH (default 0.785398 = pi/4: keeps the GB-spin increment on the principal root of the pi-periodic spin equation; never binds a physical increment ~1e-3 rad; 0 = off = old), SPINMAX (default 0 = off; >0 bounds |dw*| per step, rad), ONFAIL (retry|stop|commit, default retry), MAXITRETRY (80), CKROLL (default 20 steps; rolling pre-crash checkpoints <tag>_ckroll{0,1}.txt), EXACTEXP (1).
#   SYMFOLD (default 0 = old; 1 = KWC on the symmetry-folded orientation with a minimal-gradient lift, 2 = complex order parameter exp(i n theta)), SYMN (2: period 180 deg), GNDCURL (default 0 = GND from grad theta; 1 = from Curl Fp, Cermelli-Gurtin); tag suffix _sf<n> / _gc.
#   EOS (default vinet; murnaghan | svk = old SVK), KP (5.0), K0 (default cubic bulk modulus 137.1 GPa), VOLTAN (1: volumetric tangent correction).
#   NEWTON (default full: tangent rebuilt every global iteration + line search; modified = old), TOLRES (1e-6), RESFLOOR (1e-8), LSMAX (6), MAXIT (40). NREF only matters for NEWTON=modified.
#   TANGENT (default consistent: element tangent incl. finite strain/EOS and the implicit slip+spin response; approx = old), POLISH (20 coupled local Newton its), FDK, FDQ.
#   LINSOLVER (default gmres: GMRES + ILU(0), falls back to umfpack per solve or if the IncompleteCholesky plugin is missing), LINTOL (1e-10), PRECEVERY (0: ILU once per step), TANEVERY (1), SLIPPRED (0).
#   DISSIP (default 1: heat from the implicit dissipation tau_{n+1} dgamma >= 0; 0 = old midpoint), CDSTORE (default 1: storage capped by the plastic work; 0 = off), SOFTMIN (0.01).
# Tag pistM<DRAGRHO>v<VC>[G][R]_up<speed>ms (G: GNDPHI2=0, R: RECOV!=0), e.g. pistM3v2_up400ms, so no earlier ladder is overwritten.
# RELAXFIX=1 applies the end-of-relaxation reset (rhoS=rho0, thR=th, rhoG=0, T=T0, V=0, E0) that the relaxed checkpoint was saved before. Confined piston (cold-spray analogue); variant A resumes from out_A/A_dyn_relaxed_ck.txt, variant F starts fresh. Restart-safe, with the STALLED guard.
# usage: bash run_piston.sh <up_m_per_s> [tend_ns=30] [ramp_ns=0.5] [hold_ns=2.0] [vtkevery=10] [variant=A|F]
#   variant F (no phase field, no relaxation stage) starts from scratch; tag pistF_up<speed>ms instead of pist_up<speed>ms
#   piston velocity u_p (m/s) is applied to the right edge (u_x = -d(t)); wall u_x=0 on the left; u_y=0 on top and bottom.
#   total piston travel = u_p*(ramp+hold); frames every vtkevery*dt = 0.1 ns by default (dt = 0.01 ns).
set -uo pipefail
UPMS=${1:?piston velocity in m/s}; TEND=${2:-30}; TR=${3:-0.5}; TH=${4:-2.0}; VE=${5:-10}; V=${6:-A}
UP=$(awk -v u=$UPMS 'BEGIN{printf "%g", u/1000}')           # um/ns
ANG=${ANGLE:-45}
if [ "$ANG" = confined ]; then DIRARGS="-dir piston -pistonup $UP"; ATAG=""
else
  awk -v a="$ANG" 'BEGIN{exit !(a+0==a && a>=0 && a<=90)}' || { echo "ANGLE must be 0..90 (degrees from -x) or confined"; exit 1; }
  UN=$(awk -v u=$UP -v a=$ANG 'BEGIN{v=u*cos(a*atan2(0,-1)/180); if (v<1e-12) v=0; printf "%.10g", v}')   # normal (compression) um/ns
  UT=$(awk -v u=$UP -v a=$ANG 'BEGIN{v=u*sin(a*atan2(0,-1)/180); if (v<1e-12) v=0; printf "%.10g", v}')   # tangential (+y) um/ns
  DIRARGS="-dir shearpiston -pistonup $UN -pistonut $UT"; ATAG="_a${ANG}"
fi
if [ "$V" = F ]; then echo "mobile-drag variant: variant A only"; exit 1; fi
GP=${GNDPHI2:-1}; RC=${RECOV:-2}; CDV=${CD:-$([ "${RECOV:-2}" = 2 ] && echo 100 || echo 10)}; DR=${DRAGRHO:-3}; VCV=${VC:-2}; SFX="M${DR}v${VCV}"; [ "$GP" = 0 ] && SFX=${SFX}G; [ "$RC" != 0 ] && SFX=${SFX}R
{ [ "${DISSIP:-1}" != 0 ] || [ "${CDSTORE:-1}" != 0 ]; } && SFX=${SFX}c   # c: thermodynamically consistent heat source (implicit dissipation, storage cap)
EOSN=${EOS:-vinet}; case $EOSN in vinet) SFX=${SFX}V;; murnaghan) SFX=${SFX}Mu;; svk) ;; *) echo "EOS must be vinet|murnaghan|svk"; exit 1;; esac
DT=${DT:-0.01}; NREF=${NREF:-20}   # time step (ns) and tangent-refresh interval (steps); frames/checkpoints keep the same spacing in ns
SC=$(awk -v d=$DT 'BEGIN{printf "%d", 0.01/d + 0.5}'); [ "$SC" -ge 1 ] || SC=1; VE=$((VE*SC)); CKE=$((50*SC))
[ "$DT" != 0.01 ] || [ "$NREF" != 20 ] && SFX=${SFX}_dt${DT}_nr${NREF}   # never resume a run with a different dt
SM=${SPINMAX:-0}; [ "$SM" != 0 ] && SFX=${SFX}_sm${SM}; [ "${SYMFOLD:-0}" != 0 ] && SFX=${SFX}_sf${SYMFOLD}; [ "${GNDCURL:-0}" != 0 ] && SFX=${SFX}_gc   # sf<n>: crystal-symmetry-folded KWC orientation (-symfold), gc: GND from Curl Fp (-gndcurl)
[ -n "${RUNTAG:-}" ] && SFX=${SFX}_${RUNTAG}   # free label to keep a test series apart
[ "${SHARETAG:-0}" = 1 ] || SFX=${SFX}_mpi        # separate from the serial runs unless SHARETAG=1
TAG=pist${SFX}${ATAG}_up${UPMS}ms; D=out_$TAG; B=$TAG; mkdir -p $D
export FF_LOADPATH=${FF_LOADPATH:-$(dirname $(find $(dirname $(dirname $(which FreeFem++)))/lib* -name msh3.so 2>/dev/null | head -1) 2>/dev/null)}
done_p() { [ -s $D/${B}_log.csv ] && tail -1 $D/${B}_log.csv | cut -d, -f2 | awk -v e=$TEND -v d=$DT '{exit !($1+0==$1 && $1>=e-1.1*d)}'; }
done_p && { echo "already complete"; exit 0; }
[ -e $D/${B}_TSTOP.txt ] && { echo "ended at the temperature limit:"; cat $D/${B}_TSTOP.txt; exit 0; }
[ -e $D/${B}_NOCONV.txt ] && { echo "stopped earlier on a non-converged step:"; cat $D/${B}_NOCONV.txt; echo "(delete $D/${B}_NOCONV.txt to retry from its checkpoint, e.g. with MAXITRETRY larger or ONFAIL=commit; a non-zero SPINMAX or another DT is a new tag and starts fresh)"; exit 1; }
[ -e $D/${B}_NAN.txt ] && { echo "stopped earlier on NaN:"; cat $D/${B}_NAN.txt; echo "(last stage: $(sed -n 2p $D/${B}_stage.txt 2>/dev/null)); delete $D/${B}_NAN.txt to retry"; exit 1; }
[ -e $D/STALLED ] && { echo "STALLED earlier (see $D/STALLED); fix the cause, then delete $D/STALLED to retry"; exit 1; }
maxstep() { [ -s $D/${B}_log.csv ] && awk -F, '$1+0==$1 && $1>m {m=$1} END{print m+0}' $D/${B}_log.csv || echo 0; }
S0=$(maxstep)
if [ "$V" = A ] && [ ! -s $D/${B}_ck.txt ]; then
  [ -s out_A/A_dyn_relaxed_ck.txt ] || { echo "FAILED: out_A/A_dyn_relaxed_ck.txt missing (run variant A first)"; exit 1; }
  cp out_A/A_dyn_relaxed_ck.txt $D/${B}_ck.txt
fi
R=""; [ -s $D/${B}_ck.txt ] && R="-restart $D/${B}_ck.txt"
# resuming: drop log/diag/thermo rows written after the checkpoint (they are recomputed), so no step appears twice
if [ -n "$R" ]; then NCK=$(head -1 $D/${B}_ck.txt)
  for f in $D/${B}_log.csv $D/${B}_diag.csv $D/${B}_thermo.csv; do
    [ -s "$f" ] && awk -F, -v n=$NCK 'NR==1 && $1=="step" {print; next} $1+0<=n' "$f" > "$f.tmp" && mv "$f.tmp" "$f"
  done
  echo "resuming from the checkpoint at step $NCK (furthest step reached so far: $S0)"
fi
touch $D/$B.pvd
HDR="step,t_ns,load_GPa,iters,locres,KE_trans,KE_rot,E_elastic,E_stored,Q_heat,W_ext,D_cum,audit_residual,max_abs_slip,mean_abs_slip,max_dT,max_abs_dtheta_deg,mean_rhoS_um2,max_rhoG_um2,min_phi,neg_D_count"
[ -n "$R" ] && { [ -s $D/${B}_log.csv ] || echo "$HDR" > $D/${B}_log.csv; }   # fresh F start: the solver writes its own header
ARGS="-variant $V -j 0.01 -amp 1 $DIRARGS -ptr $TR -pth $TH -tend $TEND -dt $DT -muc 75.4 -nref $NREF -vtkevery $VE -ckevery $CKE -tstop ${TSTOP:-1400} -f0 0.005944 -kwca 0.10 -kwcs 1.05 -kwce 0.004 -kwcd 0.005 -gndphi2 $GP -recov $RC -ca ${CA:-3.16227766} -cd $CDV -diag ${DIAG:-1} -dragrho $DR -vc $VCV -rhomsat ${RHOMSAT:-100} -c2m ${C2M:-10} -rhofcut ${RHOFCUT:-1e4} -rhommin ${RHOMMIN:-1} -relaxfix ${RELAXFIX:-1} -exactexp ${EXACTEXP:-1} -spinmax $SM -spinbranch ${SPINBRANCH:-0.785398} -bckin ${BCKIN:-1} -predictor ${PREDICTOR:-vel} -onfail ${ONFAIL:-retry} -maxitretry ${MAXITRETRY:-80} -ckroll ${CKROLL:-20} -dissip ${DISSIP:-1} -cdstore ${CDSTORE:-1} -softmin ${SOFTMIN:-0.01} -eos $EOSN -kp ${KP:-5.0} -k0 ${K0:-0} -voltan ${VOLTAN:-1} -newton ${NEWTON:-full} -tolres ${TOLRES:-1e-6} -resfloor ${RESFLOOR:-1e-8} -lsmax ${LSMAX:-6} -maxit ${MAXIT:-40} -nlocal ${NLOCAL:-30} -tangent ${TANGENT:-consistent} -polish ${POLISH:-20} -fdk ${FDK:-1e-6} -fdq ${FDQ:-1e-8} -linsolver ${LINSOLVER:-gmres} -lintol ${LINTOL:-1e-10} -precevery ${PRECEVERY:-0} -tanevery ${TANEVERY:-1} -slippred ${SLIPPRED:-0} -symfold ${SYMFOLD:-0} -symn ${SYMN:-2} -gndcurl ${GNDCURL:-0}"
EDP=cosserat2d_mob_mpi.edp
[ -s $EDP ] || { echo "FAILED: $EDP missing (generate it: python3 mpi/make_mpi.py cosserat2d_mob.edp cosserat2d_mob_mpi.edp)"; exit 1; }
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
PARGS=(); [ -n "${PETSC:-}" ] && PARGS+=(-petsc "$PETSC"); [ -n "${PETSCFALL:-}" ] && PARGS+=(-petscfall "$PETSCFALL"); [ -n "${PETSCP1:-}" ] && PARGS+=(-petscp1 "$PETSCP1")
echo "[$(date +%T)] $MPIRUN FreeFem++-mpi $EDP ($NP ranks)" >> $D/$B.out
L0=$(wc -l < $D/$B.out)
$MPIRUN FreeFem++-mpi -nw -v 0 $EDP $ARGS "${PARGS[@]}" -npexpect $NP -out $D/$B $R >> $D/$B.out 2>&1
tail -n +$((L0+1)) $D/$B.out | grep -q "LAUNCHER ERROR" && { echo "launcher mismatch (see $D/$B.out): not a solver failure; set MPIRUN and resubmit"; grep -m1 "LAUNCHER ERROR" $D/$B.out; exit 1; }
RC=$?
[ "$RC" != 0 ] && [ -s $D/${B}_stage.txt ] && { echo "solver exit code $RC; last stage reached:"; cat $D/${B}_stage.txt; }
S1=$(maxstep)
if done_p; then echo COMPLETE
elif [ -e $D/${B}_TSTOP.txt ]; then echo "ENDED at the temperature limit:"; cat $D/${B}_TSTOP.txt
elif [ -e $D/${B}_NAN.txt ]; then echo "STOPPED on NaN:"; cat $D/${B}_NAN.txt
elif [ -e $D/${B}_NOCONV.txt ]; then echo "STOPPED on a non-converged step:"; cat $D/${B}_NOCONV.txt
elif [ "$S1" -le "$S0" ]; then
  { echo "STALLED: no new steps beyond step $S0 (solver exit code $RC) at $(date)"; echo "--- last lines of $D/$B.out ---"; tail -n 25 $D/$B.out; } | tee $D/STALLED
  exit 1
else echo "INCOMPLETE at step $S1 (exit code $RC): resubmit (restart-safe)"; fi
