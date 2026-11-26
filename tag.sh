# tag.sh -- sourced by every script: run name and full option list.
# in:  ANGLE (0..90, default 45), UPMS (plate speed m/s), MODE (mpi|serial), DT, RUNTAG, SHARETAG, physics/solver variables below
#      Production (version 8.0, frozen): the defaults below are the kt3pcg set-up (exact primal-dual KWC tangent -kwctan 3, Cosserat
#      multiplier -almc 4; the projected-CG phase solve and the fast assembly are solver defaults outside this key), so a run without
#      physics variables is the production run.
#      Both set-ups run the same solver (cosserat_kwc_v2*.edp) with the same numerics and parameters (Cosserat constraint by a multiplier
#      per element, mu_c 754 its augmentation; delta, tolerances, safeguards), except the micro-inertia j and therefore the theta integrator:
#      production (default): j = 0, theta quasi-static; tag suffix _v2, relaxation checkpoint relax_checkpoint_v2.txt
#      OLD=1 (or INERTIA=1): j = 0.01 um^2, theta by Newmark (average acceleration) like u; tag suffix _j, relax_checkpoint_j.txt
#      (folders ending in _old are the former OLD=1 runs with the Cosserat penalty; folders without a suffix the former base solver:
#      neither is resumed)
# out: TAG = a<ANGLE>_up<UPMS>ms<SUFFIX>, SUFFIX = [_dt<DT>]{_v2|_j}[_<RUNTAG>][_serial], UN/UT (um/ns), PHYS (the key stored in the .args files)
# A resume with other PHYS values is refused: change values together with a new RUNTAG.
ANG=${ANGLE:-45}
awk -v a="$ANG" 'BEGIN{exit !(a+0==a && a>=0 && a<=90)}' || { echo "ANGLE must be 0..90 (degrees from -x; 0 = compression, 90 = shear)" >&2; return 1 2>/dev/null || exit 1; }
DT=${DT:-0.01}
SUFFIX=""; [ "$DT" != 0.01 ] && SUFFIX="${SUFFIX}_dt${DT}"
OLD=${OLD:-${INERTIA:-0}}
case "$OLD" in 0|1) ;; *) echo "OLD must be 0 or 1" >&2; return 1 2>/dev/null || exit 1;; esac
[ "${EXP:-0}" = 1 ] && [ "$OLD" = 0 ] && echo "note: EXP=1 is no longer needed: the former experimental set-up is the production default (v2)" >&2
[ "${EXP:-0}" = 1 ] && [ "$OLD" = 1 ] && { echo "EXP=1 and OLD=1 contradict each other" >&2; return 1 2>/dev/null || exit 1; }
V2=$((1 - OLD))
KWCD=${KWCD:-5e-4}
MUC=${MUC:-754}; ALMC=${ALMC:-4}   # multiplier (mu_c = augmentation), both set-ups
if [ "$V2" = 1 ]; then SUFFIX="${SUFFIX}_v2"; J=${J:-0}      # theta quasi-static
else                   SUFFIX="${SUFFIX}_j";  J=${J:-0.01}   # rotational inertia: theta by Newmark
fi
[ -n "${RUNTAG:-}" ] && SUFFIX="${SUFFIX}_${RUNTAG}"
[ "${MODE:-mpi}" = serial ] && [ "${SHARETAG:-0}" != 1 ] && SUFFIX="${SUFFIX}_serial"
TAG=""; UN=""; UT=""
if [ -n "${UPMS:-}" ]; then
  TAG=a${ANG}_up${UPMS}ms${SUFFIX}
  UP=$(awk -v u=$UPMS 'BEGIN{printf "%g", u/1000}')                                                     # um/ns
  UN=$(awk -v u=$UP -v a=$ANG 'BEGIN{v=u*cos(a*atan2(0,-1)/180); if (v<1e-12) v=0; printf "%.10g", v}')   # compression (-x)
  UT=$(awk -v u=$UP -v a=$ANG 'BEGIN{v=u*sin(a*atan2(0,-1)/180); if (v<1e-12) v=0; printf "%.10g", v}')   # shear (+y)
fi
PHYS="-nrelax ${NRELAX:-100} -j $J -muc $MUC -f0 ${F0:-0.005944} -kwca ${KWCA:-0.10} -kwcs ${KWCS:-1.05} -kwcd $KWCD"
PHYS="$PHYS -etaphi ${ETAPHI:-0.1} -qphi ${QPHI:-100e3} -etaphimin ${ETAPHIMIN:-1e-2} -etastar ${ETASTAR:-1} -cd ${CD:-100} -ca ${CA:-3.16227766}"
PHYS="$PHYS -vc ${VC:-2} -rhomsat ${RHOMSAT:-100} -c2m ${C2M:-10} -rhofcut ${RHOFCUT:-1e4} -rhommin ${RHOMMIN:-1} -kp ${KP:-5.0} -k0 ${K0:-0} -softmin ${SOFTMIN:-0.01}"
PHYS="$PHYS -maxit ${MAXIT:-40} -maxitretry ${MAXITRETRY:-80} -tolres ${TOLRES:-1e-6} -resfloor ${RESFLOOR:-1e-8} -lsmax ${LSMAX:-6} -nlocal ${NLOCAL:-30}"
PHYS="$PHYS -polish ${POLISH:-20} -lintol ${LINTOL:-1e-10} -tanevery ${TANEVERY:-1} -precevery ${PRECEVERY:-0}"
PHYS="$PHYS -spinbranch ${SPINBRANCH:-0.785398} -jmin ${JMIN:-0.1} -stagres ${STAGRES:-100} -phitol ${PHITOL:-1e-10} -phikkt ${PHIKKT:-1e-8} -phimaxit ${PHIMAXIT:-40}"
PHYS="$PHYS -diag ${DIAG:-1} -ckroll ${CKROLL:-20}"
PHYS="$PHYS -kwctan ${KWCTAN:-3} -almc $ALMC"   # KWC tangent: 3 = exact primal-dual Newton (production), 2 = primal-dual symmetrised; Cosserat coupling: 4 = multiplier, 0 = penalty

# Kocks-Mecking constants (solver defaults K 0.05, d 5e-4 um): appended only when changed, so the default key does not move
if [ "${KMTK:-0.05}" != 0.05 ] || [ "${KMTD:-5e-4}" != 5e-4 ]; then PHYS="$PHYS -kmtk ${KMTK:-0.05} -kmtd ${KMTD:-5e-4}"; fi
