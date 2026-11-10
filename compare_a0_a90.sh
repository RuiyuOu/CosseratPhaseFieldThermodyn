#!/bin/bash
# Compare the new (Gauss + KKT phase field, -safe 1) 400 m/s runs with the earlier production run, per angle, at common times.
# usage: bash compare_a0_a90.sh [OLD_PREFIX] [NEW_PREFIX]
#   OLD_PREFIX default out_pistM3v2GRcV_gc_ne_pm1e-2_es1_mpi     (dirs  <prefix>_a<ANGLE>_up400ms/<same name without out_>_log.csv)
#   NEW_PREFIX default out_pistM3v2GRcVQ3_gc_ne_pm1e-2_es1_mpi
#   env: ANGLES ("0 90"), TIMES (ns, default "0.5 1 2 5 10 20 30")
# The new runs are NOT expected to coincide with the old ones (different element integration, exact curl, safeguards); the table shows how far apart the main quantities are
# (new/old) at each time, the last common time, and the new run's safeguard counters.  Nothing is written except to stdout.
cd "$(dirname "$0")"
OLDP=${1:-out_pistM3v2GRcV_gc_ne_pm1e-2_es1_mpi}; NEWP=${2:-out_pistM3v2GRcVQ3_gc_ne_pm1e-2_es1_mpi}
for A in ${ANGLES:-0 90}; do
  python3 - "$OLDP" "$NEWP" "$A" "${TIMES:-0.5 1 2 5 10 20 30}" <<'PY'
import sys,csv,os,glob
op,np_,A,times=sys.argv[1],sys.argv[2],sys.argv[3],[float(x) for x in sys.argv[4].split()]
def load(prefix):
    d=f"{prefix}_a{A}_up400ms"; f=os.path.join(d,d[4:]+"_log.csv")
    if not os.path.exists(f): return d,None,None
    r=list(csv.reader(open(f))); h=r[0]; rows=[]
    for l in r[1:]:
        try: rows.append(dict(zip(h,[float(x) for x in l])))
        except ValueError: pass
    return d,h,rows
do,ho,ro=load(op); dn,hn,rn=load(np_)
print(f"\n===== angle {A} deg =====\n old: {do}  {'(missing)' if ro is None else str(len(ro))+' rows, t_end '+str(ro[-1]['t_ns'])}\n new: {dn}  {'(missing)' if rn is None else str(len(rn))+' rows, t_end '+str(rn[-1]['t_ns'])}")
if not ro or not rn: sys.exit()
cols=["load_GPa","E_elastic","E_stored","Q_heat","max_dT","max_abs_slip","mean_rhoS_um2","max_rhoG_um2","min_phi","max_abs_dtheta_deg","audit_residual","iters"]
cols=[c for c in cols if c in ho and c in hn]
def at(rows,t):
    return min(rows,key=lambda r:abs(r["t_ns"]-t))
tc=min(ro[-1]["t_ns"],rn[-1]["t_ns"])
print(f" common time range: 0 .. {tc:g} ns (columns: old | new | new/old - 1)")
for t in [x for x in times if x<=tc+1e-9]+([tc] if tc not in times else []):
    a,b=at(ro,t),at(rn,t); print(f" t = {a['t_ns']:g} ns (new row t = {b['t_ns']:g})")
    for c in cols:
        x,y=a[c],b[c]; rel=(y/x-1) if abs(x)>1e-12 else float('nan')
        print(f"    {c:20s} {x:14.6g} {y:14.6g} {rel:+10.3%}")
nd=[r for r in rn if 'audit_full' in r]
if nd: print(f" new audit_full (energy balance incl. plate work), last row: {nd[-1]['audit_full']:.3e}   max |audit_full|: {max(abs(r['audit_full']) for r in nd if r['audit_full']==r['audit_full']):.3e}")
m=glob.glob(os.path.join(dn,"*_meta.json"))
if m:
    import json
    j=json.load(open(m[0])); print(" new run safeguard counters:",json.dumps(j.get("safe_mode",{})))
PY
done
