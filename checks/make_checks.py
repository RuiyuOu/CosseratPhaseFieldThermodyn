#!/usr/bin/env python3
"""Generate checks/cosserat_kwc_check.edp: the production solver plus its two self-tests (never used in production runs).

  -tancheck 1   once, at the first Newton iteration of the first step >= -tancheckstep: compare the assembled tangent with central
                differences of the global residual along a fixed direction (prints TANGENT CHECK ..., appends <out>_tancheck.txt);
  -tancheck 2   also scans the finite-difference step (the error should follow h^2 down to a round-off floor) and checks the KWC
                curvature tangent alone;
  -check 1      after the run: driving forces tau_a, tau* against -dE/d(increment) at the most stressed element, and the internal
                force against the energy derivative along a fixed direction.

Usage (from the repository root):  python3 checks/make_checks.py cosserat_kwc.edp checks/cosserat_kwc_check.edp
then see checks/run_tangent_check.sh.  Serial only.  Each insertion is an exact text match and fails loudly if the solver changed.
"""
import sys
src, dst = sys.argv[1], sys.argv[2]
s = open(src).read()
def rep(old, new, count=1):
    global s
    n = s.count(old)
    if n != count: raise SystemExit("check patch expected %d match(es), found %d:\n%s" % (count, n, old[:200]))
    s = s.replace(old, new)

TANCHECK = r'''// ---- self-test (checks/make_checks.py): tangent against central differences of the global residual ----
int tanCheck = atoi(arg("-tancheck", "0")), tanChecked = 0;
int tanCheckStep = atoi(arg("-tancheckstep", "0"));   // the check runs at the first Newton iteration of the first step >= this one (a state with slip)
func int checkTangent() {
  real[int] d(Wh.ndof); for (int i = 0; i < Wh.ndof; i++) d[i] = sin(0.37*i + 1.3); d = d .* mask;
  real eps = 1e-7*max(U.linfty, 1e-3)/max(d.linfty, 1e-30);
  real[int] U0 = U, rp(Wh.ndof), rm(Wh.ndof), dd = d; dd *= eps;
  U = U0; U += dd; u1[] = U; localUpdate(); residual(); rp = rhs;
  U = U0; U -= dd; u1[] = U; localUpdate(); residual(); rm = rhs;
  U = U0; u1[] = U; localUpdate(); residual();
  real[int] fd = rm - rp; fd *= 0.5/eps; fd = fd .* mask;          // = (d f_int/dU) d  (rhs = -f_int + ...)
  real[int] Kd = Aeff*d; Kd = Kd .* mask;
  real[int] e = Kd - fd;
  cout << "TANGENT CHECK state: max |dgamma| " << dpa.linfty << ", " << dpb.linfty << ", " << dpc.linfty << "; max |dw*| " << dwst.linfty << endl;
  cout << "TANGENT CHECK (exact Gauss-point tangent): |K d - FD| / |FD| = " << sqrt(e'*e)/max(sqrt(fd'*fd), 1e-300) << "   |FD| = " << sqrt(fd'*fd) << endl;
  { ofstream fc(out + "_tancheck.txt", append); fc << "step " << curStep << " it " << curIt << " rel.error " << sqrt(e'*e)/max(sqrt(fd'*fd), 1e-300) << endl; }
  if (tanCheck >= 2) {
    for (int ie = 0; ie < 4; ie++) {   // scan of the FD step: the error should follow h^2 down to the round-off floor
      real fac2 = 1e-5*10.^(-ie);
      real eps2 = fac2*max(U.linfty, 1e-3)/max(d.linfty, 1e-30);
      real[int] rp2(Wh.ndof), rm2(Wh.ndof), dd2 = d; dd2 *= eps2;
      U = U0; U += dd2; u1[] = U; localUpdate(); residual(); rp2 = rhs;
      U = U0; U -= dd2; u1[] = U; localUpdate(); residual(); rm2 = rhs;
      U = U0; u1[] = U; localUpdate(); residual();
      real[int] fd2 = rm2 - rp2; fd2 *= 0.5/eps2; fd2 = fd2 .* mask; real[int] e2 = Kd - fd2;
      cout << "TANGENT CHECK step scan: FD factor " << fac2 << "  |K d - FD|/|FD| = " << sqrt(e2'*e2)/max(sqrt(fd2'*fd2), 1e-300) << endl;
    }
    {   // KWC force only: FD of the couple-stress force along the theta part of d against the KWC tangent H
      Wh [dA1, dA2, dA3]; dA1[] = d; Wh [dB1, dB2, dB3]; [dB1, dB2, dB3] = [0., 0., dA3]; real[int] dth = dB1[];
      real[int] fk1(Wh.ndof), fk2(Wh.ndof); real epsk = 1e-7*max(U.linfty, 1e-3)/max(dth.linfty, 1e-30); real[int] dk = dth; dk *= epsk;
      U = U0; U += dk; u1[] = U; kin(); kwcStress(); fk1 = Bg1*mxk; rpl = Bg2*myk; fk1 += rpl;
      U = U0; U -= dk; u1[] = U; kin(); kwcStress(); fk2 = Bg1*mxk; rpl = Bg2*myk; fk2 += rpl;
      U = U0; u1[] = U; kin(); taus(); kwcStress();
      real[int] fdk = fk1 - fk2; fdk *= 0.5/epsk; fdk = fdk .* mask;
      matrix Kvc = vKv(Wh, Wh); real[int] Kdk = Kvc*dth; Kdk = Kdk .* mask; real[int] ek = Kdk - fdk;
      cout << "TANGENT CHECK KWC force only: |H d - FD|/|FD| = " << sqrt(ek'*ek)/max(sqrt(fdk'*fdk), 1e-300) << "   |FD| = " << sqrt(fdk'*fdk) << endl;
    }
  }
  return 0;
}
func int newtonSolve('''
rep('func int newtonSolve(', TANCHECK)
rep("    rn = sqrt(rhs'*rhs); if (it == 0) r0 = rn;\n", '''    if (tanCheck && !tanChecked && curStep >= tanCheckStep) { tanChecked = 1; mark(tg + "tangent check"); checkTangent(); }
    rn = sqrt(rhs'*rhs); if (it == 0) r0 = rn;
''')

ECHECK = r'''
// ---- self-test (checks/make_checks.py): -check 1 ----
if (atoi(arg("-check", "0"))) {
  // (1) tau_a = -(1/A_e) dE/dgamma_a and tau* = -(1/A_e) dE/dw* at the most stressed element (finite differences)
  localUpdate(); taus();
  real[int] q = abs(ta); int k = q.imax; real ar = Th[k].area;
  for (int sys = 0; sys < 4; sys++) {
    real h = 1e-7;
    real[int] sva = dpa, svb = dpb, svc = dpc, svw = dwst;
    if (sys == 0) dpa[k] += h; if (sys == 1) dpb[k] += h; if (sys == 2) dpc[k] += h; if (sys == 3) dwst[k] += h;
    taus(); real Ep = elasticEnergy();
    dpa = sva; dpb = svb; dpc = svc; dwst = svw;
    if (sys == 0) dpa[k] -= h; if (sys == 1) dpb[k] -= h; if (sys == 2) dpc[k] -= h; if (sys == 3) dwst[k] -= h;
    taus(); real Em = elasticEnergy();
    dpa = sva; dpb = svb; dpc = svc; dwst = svw;
    taus();
    real tt = (sys == 0) ? ta[k] : ((sys == 1) ? tb[k] : ((sys == 2) ? tcc[k] : tstar[k]));
    cout << "FD check " << (sys < 3 ? "slip " + sys : "spin  ") << ": -dE/dp/area = " << -(Ep - Em)/(2*h)/ar << "  driving force = " << tt << endl;
  }
  // (2) internal force (Gauss-point force + KWC couple stress) = dE/dU along a direction
  real[int] dirv(Wh.ndof); for (int i = 0; i < Wh.ndof; i++) dirv[i] = sin(0.37*i + 1.3)*1e-3; dirv = dirv .* mask;
  kin(); taus();
  gaussInt(1); real[int] fint = fG;
  kwcStress(); rpl = Bg1*mxk; fint += rpl; rpl = Bg2*myk; fint += rpl;
  real hU = 1e-4; real[int] U0 = U, dU = dirv; dU *= hU;
  U = U0; U += dU; u1[] = U; kin(); taus(); real Ep = elasticEnergy();
  U = U0; U -= dU; u1[] = U; kin(); taus(); real Em = elasticEnergy();
  U = U0; u1[] = U; kin(); taus();
  cout << "FD check internal force: dE/dU.d = " << (Ep - Em)/(2*hU) << "   f_int.d = " << (fint'*dirv) << endl;
}
'''
s = s.rstrip('\n') + '\n' + ECHECK
rep('// cosserat_kwc.edp  --', '// cosserat_kwc_check.edp (generated by checks/make_checks.py: production solver + -tancheck / -check self-tests)  --')
open(dst, 'w').write(s)
print("wrote", dst)
