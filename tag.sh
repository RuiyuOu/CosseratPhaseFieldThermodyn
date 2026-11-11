# tag.sh -- sourced by every script: run name and full option list.
# in:  ANGLE (0..90, default 45), UPMS (plate speed m/s), MODE (mpi|serial), DT, RUNTAG, SHARETAG, physics/solver variables below
# out: TAG = a<ANGLE>_up<UPMS>ms<SUFFIX>, SUFFIX = [_dt<DT>][_<RUNTAG>][_serial], UN/UT (um/ns), PHYS (the key stored in the .args files)
# A resume with other PHYS values is refused: change values together with a new RUNTAG.
ANG=${ANGLE:-45}
awk -v a="$ANG" 'BEGIN{exit !(a+0==a && a>=0 && a<=90)}' || { echo "ANGLE must be 0..90 (degrees from -x; 0 = compression, 90 = shear)" >&2; return 1 2>/dev/null || exit 1; }
DT=${DT:-0.01}
SUFFIX=""; [ "$DT" != 0.01 ] && SUFFIX="${SUFFIX}_dt${DT}"
[ -n "${RUNTAG:-}" ] && SUFFIX="${SUFFIX}_${RUNTAG}"
[ "${MODE:-mpi}" = serial ] && [ "${SHARETAG:-0}" != 1 ] && SUFFIX="${SUFFIX}_serial"
TAG=""; UN=""; UT=""
if [ -n "${UPMS:-}" ]; then
  TAG=a${ANG}_up${UPMS}ms${SUFFIX}
  UP=$(awk -v u=$UPMS 'BEGIN{printf "%g", u/1000}')                                                     # um/ns
  UN=$(awk -v u=$UP -v a=$ANG 'BEGIN{v=u*cos(a*atan2(0,-1)/180); if (v<1e-12) v=0; printf "%.10g", v}')   # compression (-x)
  UT=$(awk -v u=$UP -v a=$ANG 'BEGIN{v=u*sin(a*atan2(0,-1)/180); if (v<1e-12) v=0; printf "%.10g", v}')   # shear (+y)
fi
PHYS="-nrelax ${NRELAX:-100} -j ${J:-0.01} -muc ${MUC:-75.4} -f0 ${F0:-0.005944} -kwca ${KWCA:-0.10} -kwcs ${KWCS:-1.05} -kwcd ${KWCD:-0.005}"
PHYS="$PHYS -etaphi ${ETAPHI:-0.1} -qphi ${QPHI:-100e3} -etaphimin ${ETAPHIMIN:-1e-2} -etastar ${ETASTAR:-1} -cd ${CD:-100} -ca ${CA:-3.16227766}"
PHYS="$PHYS -vc ${VC:-2} -rhomsat ${RHOMSAT:-100} -c2m ${C2M:-10} -rhofcut ${RHOFCUT:-1e4} -rhommin ${RHOMMIN:-1} -kp ${KP:-5.0} -k0 ${K0:-0} -softmin ${SOFTMIN:-0.01}"
PHYS="$PHYS -maxit ${MAXIT:-40} -maxitretry ${MAXITRETRY:-80} -tolres ${TOLRES:-1e-6} -resfloor ${RESFLOOR:-1e-8} -lsmax ${LSMAX:-6} -nlocal ${NLOCAL:-30}"
PHYS="$PHYS -polish ${POLISH:-20} -fdk ${FDK:-1e-6} -fdq ${FDQ:-1e-8} -lintol ${LINTOL:-1e-10} -tanevery ${TANEVERY:-1} -precevery ${PRECEVERY:-0}"
PHYS="$PHYS -spinbranch ${SPINBRANCH:-0.785398} -jmin ${JMIN:-0.1} -stagres ${STAGRES:-100} -phitol ${PHITOL:-1e-10} -phikkt ${PHIKKT:-1e-8} -phimaxit ${PHIMAXIT:-40}"
PHYS="$PHYS -diag ${DIAG:-1} -ckroll ${CKROLL:-20}"
