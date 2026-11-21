#!/usr/bin/env python3
"""Generate cosserat_kwc_check.edp (in the repository root, next to mesh.msh): the production solver plus its three self-tests (never used in production runs).

  -tancheck 1   once, at the first Newton iteration of the first step >= -tancheckstep: compare the assembled tangent with central
                differences of the global residual along a fixed direction (prints TANGENT CHECK ..., appends <out>_tancheck.txt);
  -tancheck 2   also scans the finite-difference step (the error should follow h^2 down to a round-off floor) and checks the KWC
                curvature tangent alone;
  -check 1      after the run: (1) driving forces tau_a, tau* against -dE/d(increment) at the most stressed element (informational:
                Gauss-point energy, finite exponential update); (1b) tau = M:B against -dpsi/dh for Fp -> (I + h B) Fp on every element
                (exact, element-mean state); (2) the internal force against the energy derivative along a fixed direction.
  -tancheck with the -almc 4 multiplier: also the blocks C = d f_int/d lambda, B = d abar/dU, E = d abar/d lambda against central
                differences.
  -adcheck 1    before the time loop: the AD element-map derivatives (plugin kwc2d, elemJac) against central differences of the
                element map on every element at a perturbed state, and the plugin's element-map values against the solver's; then stop.
                Expected: |FD - AD| at the FD error (~1e-7 relative), values ~1e-13.  Steps: -fdk (1e-6, kinematics), -fdq (1e-8, increments).

Usage (from the repository root):  python3 checks/make_checks.py cosserat_kwc.edp cosserat_kwc_check.edp
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
//MULT-CHECK
  return 0;
}
func int newtonSolve('''
MULTCHECK = r'''  if (almC == 4) {   // multiplier blocks of the saddle system against central differences: C = d f_int/d lambda, B = d abar/dU, E = d abar/d lambda
    real[int] cf(15*nt + 1), mj(15*nt + 1); cf = 0;
    for (int c = 0; c < 15; c++) { real[int] col = cBm(:, c); cf(c*nt:(c+1)*nt-1) = col; }
    matrix Bm = [lmI, lmJ, cf];
    for (int c = 0; c < 15; c++) { real[int] col = cCm(:, c); cf(c*nt:(c+1)*nt-1) = col; }
    mj = mask(lmJ); cf = cf .* mj; matrix Cm = [lmJ, lmI, cf]; matrix Em = [eCon];
    real[int] dl(nt); for (int k = 0; k < nt; k++) dl[k] = cos(0.71*k + 0.2);
    real[int] L0 = lamC, U0 = U, ap(nt), am(nt), rp(Wh.ndof), rm(Wh.ndof);
    real hl = 1e-6*max(lamC.linfty, 1e-2); real[int] ddl = dl; ddl *= hl;
    lamC = L0; lamC += ddl; localUpdate(); residual(); rp = rhs; ap = aMean;
    lamC = L0; lamC -= ddl; localUpdate(); residual(); rm = rhs; am = aMean;
    lamC = L0; localUpdate(); residual();
    real[int] fdC = rm - rp; fdC *= 0.5/hl; fdC = fdC .* mask; real[int] Cd = Cm*dl; Cd = Cd .* mask; real[int] eC = Cd - fdC;
    real[int] fdE = ap - am; fdE *= 0.5/hl; real[int] Ed = Em*dl; real[int] eE = Ed - fdE;
    real[int] d(Wh.ndof); for (int i = 0; i < Wh.ndof; i++) d[i] = sin(0.37*i + 1.3); d = d .* mask;
    real eps = 1e-7*max(U.linfty, 1e-3)/max(d.linfty, 1e-30); real[int] dd = d; dd *= eps;
    U = U0; U += dd; u1[] = U; localUpdate(); residual(); ap = aMean;
    U = U0; U -= dd; u1[] = U; localUpdate(); residual(); am = aMean;
    U = U0; u1[] = U; localUpdate(); residual();
    real[int] fdB = ap - am; fdB *= 0.5/eps; real[int] Bd = Bm*d; real[int] eB = Bd - fdB;
    cout << "TANGENT CHECK multiplier C (d f_int/d lambda): |C dl - FD|/|FD| = " << sqrt(eC'*eC)/max(sqrt(fdC'*fdC), 1e-300) << "   |FD| = " << sqrt(fdC'*fdC) << endl;
    cout << "TANGENT CHECK multiplier B (d abar/dU): |B d - FD|/|FD| = " << sqrt(eB'*eB)/max(sqrt(fdB'*fdB), 1e-300) << "   |FD| = " << sqrt(fdB'*fdB) << endl;
    cout << "TANGENT CHECK multiplier E (d abar/d lambda): |E dl - FD|/|FD| = " << sqrt(eE'*eE)/max(sqrt(fdE'*fdE), 1e-300) << "   |FD| = " << sqrt(fdE'*fdE) << endl;
  }
'''
rep('func int newtonSolve(', TANCHECK)
if 'func int linSolveLM()' in s: rep('//MULT-CHECK\n', MULTCHECK)
else: rep('//MULT-CHECK\n', '')
rep("    rn = sqrt(rhs'*rhs); if (it == 0) r0 = rn;\n", '''    if (tanCheck && !tanChecked && curStep >= tanCheckStep) { tanChecked = 1; mark(tg + "tangent check"); checkTangent(); }
    rn = sqrt(rhs'*rhs); if (it == 0) r0 = rn;
''')

ADFD = r'''// ---- self-test (checks/make_checks.py): central-difference element map for -adcheck ----
real fdK = atof(arg("-fdk", "1e-6")), fdQ = atof(arg("-fdq", "1e-8"));
// perturb input j of the element map: 0 F11, 1 F12 (= du1/dy), 2 F21 (= du2/dx), 3 F22, 4 theta, 5-7 dgamma_1..3, 8 dw*
func int addIn(int j, real h) {
  if (j == 0) e11 += h; if (j == 1) hy1 += h; if (j == 2) hx2 += h; if (j == 3) e22 += h; if (j == 4) thc += h;
  if (j == 5) dpa += h; if (j == 6) dpb += h; if (j == 7) dpc += h; if (j == 8) dwst += h;
  return 0;
}
// output i of the element map: 0-4 finite-strain stress (P11, P12, P21, P22, dpsi/dtheta), 5-8 tau_1..3, tau*
func int outI(int i, real[int] & o) {
  if (i == 0) o = PF11; if (i == 1) o = PF12; if (i == 2) o = PF21; if (i == 3) o = PF22; if (i == 4) o = psth;
  if (i == 5) o = ta; if (i == 6) o = tb; if (i == 7) o = tcc; if (i == 8) o = tstar;
  return 0;
}
// central differences of the element map w.r.t. inputs j0..j1 (-fdk for the kinematics, -fdq for the increments)
func int fdMapFD(int j0, int j1) {
  real[int] op(nt), om(nt), d(nt);
  for (int j = j0; j <= j1; j++) {
    real h = (j < 5) ? fdK : fdQ;
    real[int,int] Sp(nt, 9), Sm(nt, 9);
    addIn(j, h); taus(); for (int i = 0; i < 9; i++) { outI(i, op); Sp(:, i) = op; }
    addIn(j, -2*h); taus(); for (int i = 0; i < 9; i++) { outI(i, om); Sm(:, i) = om; }
    addIn(j, h);
    for (int i = 0; i < 9; i++) {
      op = Sp(:, i); om = Sm(:, i); d = op - om; d *= 0.5/h;
      if (i < 5) DSF(:, i*9 + j) = d;
      else DTQ(:, (i-5)*9 + j) = d;
    }
  }
  taus();
  return 0;
}
// exact local Jacobian dr/dq'''
rep('// exact local Jacobian dr/dq', ADFD)

ADCHECK = r"""if (atoi(arg("-adcheck", "0"))) {   // self-test (checks/make_checks.py): AD vs central differences of the element map, all elements
  randinit(5);
  for (int k = 0; k < nt; k++) { e11[k] += 0.01*(2*randreal1() - 1); hy1[k] += 0.01*(2*randreal1() - 1); hx2[k] += 0.01*(2*randreal1() - 1); e22[k] += 0.01*(2*randreal1() - 1);
    thc[k] += 0.05*(2*randreal1() - 1); dpa[k] = 1e-3*(2*randreal1() - 1); dpb[k] = 1e-3*(2*randreal1() - 1); dpc[k] = 1e-3*(2*randreal1() - 1); dwst[k] = 1e-3*(2*randreal1() - 1); }
  real[int] s0 = e11, s1 = hy1, s2 = hx2, s3 = e22, s4 = thc;
  real[int,int] FS(nt, 45), FT(nt, 36); real tf0 = clock(); fdMapFD(0, 8); real tfd = clock() - tf0; FS = DSF; FT = DTQ;
  e11 = s0; hy1 = s1; hx2 = s2; e22 = s3; thc = s4;
  tf0 = clock(); elemJac(); real tad = clock() - tf0;
  real[int,int] VO(nt, 10); k2dMap(ADin, VO);
  real ev = 0; { real[int] o(nt); outI(0, o); real[int] d = VO(:, 0); d -= o; ev = max(ev, d.linfty); outI(4, o); d = VO(:, 4); d -= o; ev = max(ev, d.linfty); outI(5, o); d = VO(:, 5); d -= o; ev = max(ev, d.linfty); outI(8, o); d = VO(:, 8); d -= o; ev = max(ev, d.linfty); }
  for (int i = 0; i < 9; i++) { real e = 0, m = 0; for (int j = 0; j < 9; j++) { real[int] a(nt), f(nt);
      if (i < 5) { a = DSF(:, i*9 + j); f = FS(:, i*9 + j); } else { a = DTQ(:, (i-5)*9 + j); f = FT(:, (i-5)*9 + j); }
      m = max(m, a.linfty); f -= a; e = max(e, f.linfty); }
    cout << "adcheck: output " << i << "  max |FD - AD| " << e << "  (max |AD| " << m << ")" << endl; }
  cout << "adcheck: values plugin vs solver max |diff| " << ev << ";  time FD " << tfd << " s, AD " << tad << " s (" << nt << " elements)" << endl;
  exit(0);
}
refreshTangent();

int nsteps"""
rep("refreshTangent();\n\nint nsteps", ADCHECK)

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
    cout << "FD check " << (sys < 3 ? "slip " + sys : "spin  ") << " (finite increment, Gauss-point energy; informational): -dE/dp/area = " << -(Ep - Em)/(2*h)/ar << "  driving force = " << tt << endl;
  }
  // (1b) driving forces tau = M : B as the derivative of the energy for an infinitesimal left-multiplicative change of the trial
  // plastic distortion, Fp -> (I + h B) Fp, at the element-mean state the local update uses: tau = -dpsi/dh exactly, every element.
  // (1) differs from tau through the Gauss-point sampling of the energy and through the finite exponential update
  // (dFp/dq Fp^-1 = int_0^1 exp(sA) B exp(-sA) ds, which equals B only for A = 0 or [A, B] = 0).
  {
    kin(); taus();
    real[int] T0 = ta, T1 = tb, T2 = tcc, T3 = tstar, s11 = tf11, s12 = tf12, s21 = tf21, s22 = tf22, ps(nt), ms(nt), b11(nt), b12(nt), b21(nt), b22(nt), w(nt);
    real hb = 1e-5;   // the error is round-off (falls like 1/h: 2e-7, 2e-8, 2e-9 at 1e-7, 1e-6, 1e-5 on the coarse test)
    for (int g = 0; g < 4; g++) {
      if (g == 0) { b11 = A11a; b12 = A12a; b21 = A21a; b22 = A22a; }
      if (g == 1) { b11 = A11b; b12 = A12b; b21 = A21b; b22 = A22b; }
      if (g == 2) { b11 = A11c; b12 = A12c; b21 = A21c; b22 = A22c; }
      if (g == 3) { b11 = 0; b12 = -1; b21 = 1; b22 = 0; }   // plastic spin generator (trialFp: L12 -= dw*, L21 += dw*)
      for (int sg = 0; sg < 2; sg++) {
        real hh = (sg == 0) ? hb : -hb;
        tf11 = b11 .* s11; w = b12 .* s21; tf11 += w; tf11 *= hh; tf11 += s11;
        tf12 = b11 .* s12; w = b12 .* s22; tf12 += w; tf12 *= hh; tf12 += s12;
        tf21 = b21 .* s11; w = b22 .* s21; tf21 += w; tf21 *= hh; tf21 += s21;
        tf22 = b21 .* s12; w = b22 .* s22; tf22 += w; tf22 *= hh; tf22 += s22;
        finiteStress(); if (sg == 0) ps = psiF; else ms = psiF;
      }
      real[int] dr = ps - ms; dr *= -0.5/hb;
      real[int] tt = T0; if (g == 1) tt = T1; if (g == 2) tt = T2; if (g == 3) tt = T3;
      real[int] e = dr - tt;
      cout << "FD check Mandel " << (g < 3 ? "slip " + g : "spin  ") << " (tau = M:B, all elements): max |-dpsi/dh - tau| / max |tau| = " << e.linfty/max(tt.linfty, 1e-300) << "   max |tau| = " << tt.linfty << endl;
    }
    tf11 = s11; tf12 = s12; tf21 = s21; tf22 = s22; taus();
  }
  // (2) internal force (Gauss-point force + KWC couple stress) = dE/dU along a direction
  real[int] dirv(Wh.ndof); for (int i = 0; i < Wh.ndof; i++) dirv[i] = sin(0.37*i + 1.3)*1e-3; dirv = dirv .* mask;
  kin(); taus();
  updatePhi2(); gaussInt(1); real[int] fint = fG;   // element mean of g(phi) of the current phi: the check runs after the phase step
  kwcStress(); rpl = Bg1*mxk; fint += rpl; rpl = Bg2*myk; fint += rpl;
  real hU = 1e-4; real[int] U0 = U, dU = dirv; dU *= hU;
  U = U0; U += dU; u1[] = U; kin(); taus(); real Ep = elasticEnergy();
  U = U0; U -= dU; u1[] = U; kin(); taus(); real Em = elasticEnergy();
  U = U0; u1[] = U; kin(); taus();
  cout << "FD check internal force: dE/dU.d = " << (Ep - Em)/(2*hU) << "   f_int.d = " << (fint'*dirv) << "   relative difference = " << abs((Ep - Em)/(2*hU) - fint'*dirv)/max(abs(fint'*dirv), 1e-300) << endl;
}
'''
s = s.rstrip('\n') + '\n' + ECHECK
rep('// cosserat_kwc.edp  --', '// cosserat_kwc_check.edp (generated by checks/make_checks.py: production solver + -tancheck / -check / -adcheck self-tests)  --')
open(dst, 'w').write(s)
print("wrote", dst)
