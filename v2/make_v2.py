#!/usr/bin/env python3
"""Generate the production solver cosserat_kwc_v2.edp from cosserat_kwc.edp (exact-text patches; the base solver is not touched).
With -kwctan 0 -almc 0 -asm 0 (the defaults of the file) the v2 solver takes the same steps as the base solver (the OLD=1 set-up);
run.sh passes -kwctan 2 -almc 4 -asm 1 for production runs.

  -kwctan 2   KWC term c_e N(g_e), g_e = grad theta_t on element e: primal-dual Newton tangent (Chan, Golub, Mulet 1999).  A dual flux
              w_e (|w_e| <= 1) per element, H = (c/n)(I - sym(w g^T)/n), dual step shortened to keep |w| <= 1.  The residual and the
              convergence test are unchanged, so a converged step is the same discrete solution as with the exact tangent (-kwctan 0).
  -kwctan 3   as 2 without the symmetrisation: H = (c/n)(I - w g^T/n).  With the dual correction eliminated element by element this
              is the exact reduced Newton tangent of the primal-dual system (flux c w, dual equation n w - g = 0) and its undamped
              direction at every iteration; primal backtracking and the separately safeguarded dual update (accepted dg, |w| <= 1)
              make the update a plain Newton step only where they are inactive.  The theta-theta block is unsymmetric.
  -almc 4     Cosserat constraint abar_e = 0 (Gauss mean of a = skew(Ubar) part - theta mismatch, one multiplier lambda_e per element),
              psi = lambda a + 2 mu_c a^2, solved monolithically: lambda is a Newton unknown, each iteration solves
              [[K, C], [B, E]] (dU, dlambda) = (R, -abar).  B = d abar/dU, C = d f_int/d lambda, E = d abar/d lambda, all including the
              implicit slip/spin response (derivatives from the plugin, k2dJacL).  Serial: sparse direct solve.  MPI (mpi/make_mpi.py):
              FGMRES with a Schur preconditioner, MUMPS fallback (-almlin, -almpetsc, -almfall).  Merit function
              sqrt(|R|^2 + sum_e (s A_e abar_e)^2), s = 4 max(mu_c, 75.4); the line search moves U and lambda together.
              mu_c (-muc) is the augmentation and must be > 0.
  -asm 1      tangent matrix assembled directly by the plugin (k2dAsm*): pattern, kinematic operators, element -> entry map and constant
              part built once; per iteration only the Gauss-point blocks HQ, CQ and the KWC coefficients.  Same matrix as the forms
              vJg + vKv + Kconst; z'Ax is compared with the FreeFEM matrix after every assembly.  (run.sh: FASTASM=1)
  -asmcheck N with -asm 1: every N-th Newton iteration of the run (the first included) the matrix is also assembled by the forms and
              compared entry by entry; stops if they differ by more than 1e-12 of the largest entry (run.sh: ASMCHECK=N)
  -betamax b  lattice-rotation branch guard (radians): beta = atan2(Ub21 - Ub12, Ub11 + Ub22) per element; a trial iterate with |beta| > 2b
              is inadmissible (line-search backtracking), an accepted step needs |beta| <= b; largest |beta| per step in <out>_newton.csv
              (run.sh: pi/4)
  -slipexact 1 slip rate law 1/|gd| = 1/(g0 x^m) + eta/|tau| evaluated exactly for every x (no cap at x = 300), with the derivative of that
              expression (run.sh: on for v2)
  -localtol r  after the global Newton converges, the local update is redone and the four local residuals |r_a|, |r_*| must be below r
              (increment units); otherwise the step goes to the retry path and, failing again, stops cleanly (run.sh: 1e-10 for v2, LOCALTOL)
  -lin direct sparse direct solve of the Newton correction (serial only; checks)
  -lam0 v     start value of the multiplier (self-test of the plugin derivatives with -adcheck)
  checkpoints lambda and w are in the side file <checkpoint>.mult
  output      <out>_newton.csv: per step (it = -1 rows) iterations, wall seconds, max |abar|;
              with -nlogit 1 also one row per global iteration (step, it, merit, clock, max |abar|)

usage: python3 v2/make_v2.py cosserat_kwc.edp cosserat_kwc_v2.edp"""
import sys
src, dst = sys.argv[1:3]
s = open(src).read()
def rep(old, new, count=1):
    global s
    n = s.count(old)
    if n != count: raise SystemExit("make_v2: expected %d match(es), found %d:\n%s" % (count, n, old[:200]))
    # a // comment in the last line of the replacement must not swallow code that follows the match on the same line
    last = new.split('\n')[-1]
    if '//' in last and '"' not in last.split('//')[-1]:
        i = 0
        while True:
            i = s.find(old, i)
            if i < 0: break
            tail = s[i + len(old):].split('\n', 1)[0]
            head = s[:i].rsplit('\n', 1)[-1] + old.split('\n')[-1]   # the tail is already a comment in the base solver
            if tail.strip() and not tail.lstrip().startswith('//') and '//' not in head:
                raise SystemExit("make_v2: the replacement ends in a // comment but code follows on the same line (%r):\n%s" % (tail[:80], old[:200]))
            i += len(old)
    s = s.replace(old, new)

rep("// Model, parameters and numerics: README.md; code map: CODE_MAP.md.",
    "// PRODUCTION solver v2: generated by v2/make_v2.py from cosserat_kwc.edp (never edit; options: v2/README.md).\n"
    "// Model, parameters and numerics: README.md; code map: CODE_MAP.md.")

# ---- declarations ----
rep("Ph H11 = 0, H12 = 0, H22 = 0;   // KWC curvature tangent (refreshTangent)\n",
    r'''Ph H11 = 0, H12 = 0, H22 = 0;   // KWC curvature tangent (refreshTangent)
Ph Hd = 0;   // its antisymmetric part (H12 - H21)/2, row = flux component (-kwctan 3; 0 otherwise)
// ==== v2 options (v2/make_v2.py) ====
int kwcTan = atoi(arg("-kwctan", "0")), almC = atoi(arg("-almc", "0")), nlogIt = atoi(arg("-nlogit", "0"));
int asmFast = atoi(arg("-asm", "1"));
int asmCheck = atoi(arg("-asmcheck", "0"));   // with -asm 1: every N-th Newton iteration of the run (counted over all steps, the first one included)
                                              // the matrix is also assembled by the forms and compared entry by entry (ASMCHECK=N)
int nAsmIt = 0, nAsmChk = 0; bool asmChkNow = false, jgFresh = false; real asmChkMax = 0;
if (asmCheck < 0 || (asmCheck > 0 && !asmFast)) { cout << "ERROR: -asmcheck N needs -asm 1 (it compares the direct assembly with the forms)" << endl; exit(1); }   // 1: tangent matrix assembled directly by the plugin on a fixed pattern (same matrix as the varf assembly)
real tLs = 0, tStep0 = 0, tNewt0 = 0, tNewt1 = 0;   // whole-step timing
func real wallNow() { return clock(); }   // wall seconds of <out>_newton.csv (mpi/make_mpi.py: mpiWtime())
real asmDev = 0; bool asmReady = false, asmJok = false; real[int] asmVals(1), asmX(1), asmZ(1); real[int,int] asmHK(nt, 4);
// serial: once a direct solver has been attached to Aeff (the GMRES fallback of linSolve, or a direct linear solver), Aeff.coef no
// longer reaches the matrix the product uses, so the matrix is rebuilt from its triplets (same entry order, checked by z'Ax below)
bool asmSerialDirect = true; int[int] asmI(1), asmJ(1); int asmFallSeen = 0;
real[int,int] geoA(nt, 55); int[int,int] connA(nt, 15);   // gbX, gbY, P1 gradients, area; element dofs (filled with the operators below)
real almTol = atof(arg("-almtol", "1e-8")), almFin = atof(arg("-almfinal", "10"));   // Newton needs max |abar| <= almTol; the final state is accepted up to almFin x almTol
if (kwcTan != 0 && kwcTan != 2 && kwcTan != 3) { cout << "ERROR: -kwctan must be 0 (exact tangent), 2 (primal-dual, symmetrised) or 3 (primal-dual, exact Newton of the primal-dual system)" << endl; exit(1); }
if (almC != 0 && almC != 4) { cout << "ERROR: -almc must be 0 (penalty) or 4 (monolithic multiplier)" << endl; exit(1); }
if (almC == 4 && muc <= 0) { cout << "ERROR: -almc 4 needs -muc > 0: without the augmentation 2 mu_c a^2 the lattice rotation inside an element is not controlled and the iteration runs away" << endl; exit(1); }
real[int] wdx(nt), wdy(nt), gox(nt), goy(nt), wsx(nt), wsy(nt); bool wInit = false;            // primal-dual: dual flux w, gradient before the step, step-start copy
real[int] lamC(nt), aMean(nt), lam0(nt);                                                       // Cosserat multiplier, Gauss mean of a, step-start copy
lamC = atof(arg("-lam0", "0")); aMean = 0;   // -lam0: start value of the multiplier (self-tests)
// -almc 4 (monolithic multiplier): plugin derivatives DAL (da/d(k,q), d(stress)/dlambda, d(tau)/dlambda), dq/dlambda, coefficients of the
// constraint blocks B (rows: elements) and C (columns: elements) on the 15 element dofs, E (diagonal), Newton correction of lambda
real[int,int] DAL(nt, 18), DQL(nt, 4), cBm(nt, 15), cCm(nt, 15), accB(nt, 5); real[int] eCon(nt), dLam(nt), lamBase(nt); eCon = 0; dLam = 0;
real almScale = 4*max(muc, 75.4);   // GPa: weight of the constraint residual in the merit function
bool multRead = false;
real aMax = 0;
string nlogFile = isRoot() ? (arg("-out", "run") + "_newton.csv") : "/dev/null";
bool nlHdr = (exec("test -s " + nlogFile) != 0);   // header whenever the file is new (also for a run that starts from the relaxation checkpoint)
ofstream nlog(nlogFile, append); if (nlHdr) nlog << "step,it,res_or_iters,clock_or_stepwall,max_abs_a,max_abs_beta_deg" << endl;
''')

# ---- serial direct solve (the MPI generator overrides linMode) ----
rep('int linMode = (linName == "gmres") ? 1 : 0;\n',
    'int linMode = (linName == "gmres") ? 1 : 0;\nif (arg("-lin", "gmres") == "direct") linMode = 0;   // v2: sparse direct solve (serial)\n')

# ---- Cosserat multiplier in the element stress (solver side; the plugin gets lambda in elemJac) ----
rep("  real[int] ska = aa; ska *= 2*muc;\n",
    "  real[int] ska = aa; ska *= 2*muc;\n"
    "  if (almC) { tmp = lamC .* aa; psiF += tmp; tmp = lamC; tmp *= 0.5; ska += tmp; }   // + lambda a\n")
rep("  k2dJac(ADin, DSF, DTQ);\n", "  if (almC) k2dLam(lamC);\n  if (almC == 4) k2dJacL(ADin, DSF, DTQ, DAL); else k2dJac(ADin, DSF, DTQ);\n")

# ---- -almc 4: basis values at the Gauss points / centroid and the triplet pattern of B, C ----
rep('''    GX1[q] = GD1[q]'; GX2[q] = GD2[q]'; GX3[q] = GD3[q]'; GX4[q] = GD4[q]'; GX5[q] = GD5[q]';
  }
}
''', '''    GX1[q] = GD1[q]'; GX2[q] = GD2[q]'; GX3[q] = GD3[q]'; GX4[q] = GD4[q]'; GX5[q] = GD5[q]';
  }
}
// v2 (-almc 4): P2 basis gradients on every element at the 3 Gauss points (q = 0..2) and at the centroid (q = 3), column q*6 + i for the
// local dof i (same formulas as GD above); P1 values of theta there; (element, dof) pattern of the constraint blocks
real[int,int] gbX(nt, 24), gbY(nt, 24), gbT(4, 3); int[int] lmI(15*nt + 1), lmJ(15*nt + 1);
{
  gbT = 0; gbT(0,0) = 0.5; gbT(0,1) = 0.5; gbT(1,1) = 0.5; gbT(1,2) = 0.5; gbT(2,0) = 0.5; gbT(2,2) = 0.5; gbT(3,0) = 1./3; gbT(3,1) = 1./3; gbT(3,2) = 1./3;
  for (int k = 0; k < nt; k++) {
    real x0 = Th[k][0].x, y0 = Th[k][0].y, x1 = Th[k][1].x, y1 = Th[k][1].y, x2 = Th[k][2].x, y2 = Th[k][2].y;
    real A2 = (x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0);
    real[int] gx(3), gy(3);
    gx[0] = (y1 - y2)/A2; gy[0] = (x2 - x1)/A2; gx[1] = (y2 - y0)/A2; gy[1] = (x0 - x2)/A2; gx[2] = (y0 - y1)/A2; gy[2] = (x1 - x0)/A2;
    for (int q = 0; q < 4; q++) {
      for (int i = 0; i < 3; i++) { gbX(k, q*6 + i) = (4*gbT(q,i) - 1)*gx[i]; gbY(k, q*6 + i) = (4*gbT(q,i) - 1)*gy[i]; }
      for (int m = 0; m < 3; m++) { int a = (m + 1)%3, b = (m + 2)%3; gbX(k, q*6 + 3 + m) = 4*(gbT(q,b)*gx[a] + gbT(q,a)*gx[b]); gbY(k, q*6 + 3 + m) = 4*(gbT(q,b)*gy[a] + gbT(q,a)*gy[b]); }
    }
    for (int c = 0; c < 15; c++) { lmI[c*nt + k] = k; lmJ[c*nt + k] = Wh(k, c); connA(k, c) = Wh(k, c); }
    for (int c = 0; c < 24; c++) { geoA(k, c) = gbX(k, c); geoA(k, 24 + c) = gbY(k, c); }
    for (int i = 0; i < 3; i++) { geoA(k, 48 + i) = gx[i]; geoA(k, 51 + i) = gy[i]; }
    geoA(k, 54) = 0.5*A2;
  }
  lmI[15*nt] = nt - 1; lmJ[15*nt] = Wh.ndof - 1;   // zero pad entry fixes the matrix shape
}
''')

# ---- -almc 4: dq/dlambda (sixth right-hand side of the local solve) and the constraint blocks ----
rep('''  int mm = 5;   // right-hand sides: the 5 kinematic inputs
''', '''  int mm = (almC == 4) ? 6 : 5;   // right-hand sides: the 5 kinematic inputs (+ lambda)
''')
rep('''    real[int] col(nt); col = DTQ(:, a*9 + j);
    if (a == 0) { col = col .* dga; col *= dt; }''', '''    real[int] col(nt); if (j < 5) col = DTQ(:, a*9 + j); else col = DAL(:, 14 + a);
    if (a == 0) { col = col .* dga; col *= dt; }''')
rep('''  solve4(mm);   // BB = dq/dk (4 x 5)
''', '''  solve4(mm);   // BB = dq/dk (4 x 5) (+ dq/dlambda)
  if (almC == 4) for (int b = 0; b < 4; b++) DQL(:, b) = BB(:, b*6 + 5);
''')
rep('''    elemJac();
    for (int i = 0; i < 5; i++) for (int j = 0; j < 5; j++) {
      d = DSF(:, i*9 + j); HQ(:, q*25 + i*5 + j) = d;
''', '''    elemJac();
    if (almC == 4) {   // constraint blocks: coefficients on the element dofs, summed over the Gauss points
      if (q == 0) { cBm = 0; cCm = 0; accB = 0; eCon = 0; }
      real[int] cc(nt), gb(nt);
      for (int j = 0; j < 5; j++) {
        d = DAL(:, j); d *= gqw(q);                                     // d a_q/dk_q at fixed increments -> B through the Gauss-point kinematics
        int c0 = (j < 2) ? 0 : ((j < 4) ? 6 : 12), nc = (j < 4) ? 6 : 3;
        for (int i = 0; i < nc; i++) {
          if (j == 4) { gb = d; gb *= gbT(q, i); } else { if (j == 0 || j == 2) gb = gbX(:, q*6 + i); else gb = gbY(:, q*6 + i); gb = gb .* d; }
          cc = cBm(:, c0 + i); cc += gb; cBm(:, c0 + i) = cc;
        }
        d = DAL(:, 5); t2 = BQ(:, j); d = d .* t2;                      // through the increments: d a_q/dq . dq/dkbar (centroid kinematics)
        for (int b = 1; b < 4; b++) { t1 = DAL(:, 5 + b); t2 = BQ(:, b*5 + j); t1 = t1 .* t2; d += t1; }
        d *= gqw(q); cc = accB(:, j); cc += d; accB(:, j) = cc;
      }
      for (int b = 0; b < 4; b++) { t1 = DAL(:, 5 + b); t2 = DQL(:, b); t1 = t1 .* t2; t1 *= gqw(q); eCon += t1; }   // E = d abar/d lambda
      for (int i = 0; i < 5; i++) {                                     // C = d f_int/d lambda = sum_q GX_q' A w_q (dS_q/dlambda + dS_q/dq dq/dlambda)
        d = DAL(:, 9 + i);
        for (int b = 0; b < 4; b++) { t1 = DSF(:, i*9 + 5 + b); t2 = DQL(:, b); t1 = t1 .* t2; d += t1; }
        d = d .* areas; d *= gqw(q);
        int c0 = (i < 2) ? 0 : ((i < 4) ? 6 : 12), nc = (i < 4) ? 6 : 3;
        for (int l = 0; l < nc; l++) {
          if (i == 4) { gb = d; gb *= gbT(q, l); } else { if (i == 0 || i == 2) gb = gbX(:, q*6 + l); else gb = gbY(:, q*6 + l); gb = gb .* d; }
          cc = cCm(:, c0 + l); cc += gb; cCm(:, c0 + l) = cc;
        }
      }
      if (q == 2) for (int j = 0; j < 5; j++) {                         // centroid part of B
        d = accB(:, j);
        int c0 = (j < 2) ? 0 : ((j < 4) ? 6 : 12), nc = (j < 4) ? 6 : 3;
        for (int i = 0; i < nc; i++) {
          if (j == 4) { gb = d; gb *= gbT(3, i); } else { if (j == 0 || j == 2) gb = gbX(:, 18 + i); else gb = gbY(:, 18 + i); gb = gb .* d; }
          cc = cBm(:, c0 + i); cc += gb; cBm(:, c0 + i) = cc;
        }
      }
    }
    for (int i = 0; i < 5; i++) for (int j = 0; j < 5; j++) {
      d = DSF(:, i*9 + j); HQ(:, q*25 + i*5 + j) = d;
''')

# ---- Gauss mean of a ----
rep('''    finiteStress();
    z = psiF .* areas; gaussE += gqw(q)*z.sum;
''', '''    finiteStress();
    if (almC) { if (q == 0) aMean = 0; tmp = aa; tmp *= gqw(q); aMean += tmp; }
    z = psiF .* areas; gaussE += gqw(q)*z.sum;
''')

# ---- KWC tangent ----
rep('''    real[int] nr = gx .* gx; tmp = gy .* gy; nr += tmp; nr += dlt*dlt; nr = sqrt(nr);
    // consistent KWC tangent: f0*s*phi2/n*(I - g g/n^2) + f0*s*phi2*dlt^2/n^3 * g g/n^2 ... written as
    // H = a I + c g(x)g,  a = f0(s gbar/n + e^2),  c = -f0 s gbar/n^3   (gbar = element mean of g(phi))
    real[int] aH = gm ./ nr; aH *= sK; aH *= f0;
    real[int] n3 = nr .* nr; n3 = n3 .* nr; real[int] cc = gm ./ n3; cc *= -f0*sK;
    tmp = gx .* gx; tmp = tmp .* cc; tmp += aH; H11[] = tmp;
    tmp = gy .* gy; tmp = tmp .* cc; tmp += aH; H22[] = tmp;
    tmp = gx .* gy; tmp = tmp .* cc; H12[] = tmp;
''', '''    real[int] nr = gx .* gx; tmp = gy .* gy; nr += tmp; nr += dlt*dlt; nr = sqrt(nr);
    // consistent KWC tangent: H = a I + c g(x)g,  a = f0 s gbar/n,  c = -f0 s gbar/n^3   (gbar = element mean of g(phi))
    real[int] aH = gm ./ nr; aH *= sK; aH *= f0;
    if (kwcTan >= 2) {                                                                    // primal-dual: (a/n)(n I - sym(w g^T)) [+ skew part, 3]
      real[int] cc = aH ./ nr; cc *= -1;
      tmp = wdx .* gx; tmp = tmp .* cc; tmp += aH; H11[] = tmp;
      tmp = wdy .* gy; tmp = tmp .* cc; tmp += aH; H22[] = tmp;
      tmp = wdx .* gy; tmp2 = wdy .* gx; tmp += tmp2; tmp *= 0.5; tmp = tmp .* cc; H12[] = tmp;
      if (kwcTan == 3) { tmp = wdx .* gy; tmp2 = wdy .* gx; tmp -= tmp2; tmp *= 0.5; tmp = tmp .* cc; Hd[] = tmp; }   // (H12 - H21)/2 of (a/n)(I - w g^T/n)
    } else {
    real[int] n3 = nr .* nr; n3 = n3 .* nr; real[int] cc = gm ./ n3; cc *= -f0*sK;
    tmp = gx .* gx; tmp = tmp .* cc; tmp += aH; H11[] = tmp;
    tmp = gy .* gy; tmp = tmp .* cc; tmp += aH; H22[] = tmp;
    tmp = gx .* gy; tmp = tmp .* cc; H12[] = tmp;
    }
''')

# ---- KWC tangent form with the antisymmetric part (-kwctan 3): H_ij d_j(theta) d_i(zeta), Hd = (H12 - H21)/2 ----
rep("  int2d(Th)( H11*dx(th)*dx(ze) + H12*(dx(th)*dy(ze) + dy(th)*dx(ze)) + H22*dy(th)*dy(ze) );\n",
    "  int2d(Th)( H11*dx(th)*dx(ze) + H12*(dx(th)*dy(ze) + dy(th)*dx(ze)) + H22*dy(th)*dy(ze) + Hd*(dy(th)*dx(ze) - dx(th)*dy(ze)) );\n")

# ---- monolithic multiplier: constraint residual and saddle-point solve ----
rep("func int newtonSolve(int mx, int itOff, string tg) {\n", r"""// ==== v2: monolithic Cosserat multiplier (-almc 4) ====
func real conRes() { if (almC != 4) return 0.; real[int] qc = aMean .* areas; qc *= almScale; qc = qc .* qc; return gsumEs(qc); }
func int linSolveLM() {   // saddle-point Newton correction for (U, lambda)
  real[int] cf(15*nt + 1), mj(15*nt + 1); cf = 0;
  for (int c = 0; c < 15; c++) { real[int] col = cBm(:, c); cf(c*nt:(c+1)*nt-1) = col; }
  matrix Bm = [lmI, lmJ, cf];
  for (int c = 0; c < 15; c++) { real[int] col = cCm(:, c); cf(c*nt:(c+1)*nt-1) = col; }
  mj = mask(lmJ); cf = cf .* mj;                    // no multiplier force on the prescribed dofs
  matrix Cm = [lmJ, lmI, cf];
  matrix Em = [eCon];
  real[int] ca = aMean; ca *= -1;
  //LM-SOLVE-BEGIN
  matrix MM = [[Aeff, Cm], [Bm, Em]];
  set(MM, solver = sparsesolver);
  int nW = Wh.ndof; real[int] rr(nW + nt), xx(nW + nt);
  rr(0:nW-1) = rhs; rr(nW:nW+nt-1) = ca;
  xx = MM^-1*rr;
  tv = xx(0:nW-1); dLam = xx(nW:nW+nt-1);
  //LM-SOLVE-END
  return 0;
}
func int newtonSolve(int mx, int itOff, string tg) {
""")

# ---- global Newton: iteration counter, multiplier convergence, updates after the accepted step ----
rep("    rn = sqrt(rhs'*rhs); if (it == 0) r0 = rn;\n",
    '''    rn = sqrt(rhs'*rhs); if (it == 0) r0 = rn;
    if (almC == 4) { rn = sqrt(rn*rn + conRes()); if (it == 0) r0 = rn; }   // merit: residual + weighted constraint
    if (almC) { real[int] qa = abs(aMean); aMax = gmaxS(qa.max); }
    if (nlogIt) nlog << curStep << "," << it << "," << rn << "," << clock() << "," << aMax << "," << betaIt*180/pi << endl;
''')
rep('''      if (rn <= tolR*max(r0, rFloor)) break;                    // converged: increment and residual (local state is current)
      if (du < 0.01*tolU) {''',
    '''      bool mulOK = (!almC || aMax <= almTol);   // constraint converged too
      if (rn <= tolR*max(r0, rFloor) && mulOK) break;           // converged: increment and residual (local state is current)
      if (du < 0.01*tolU && mulOK) {''')
rep('''    mark(tg + "linear solve for (u, theta)");
    linSolve();
''', '''    mark(tg + "linear solve for (u, theta)");
    gox = gsx; goy = gsy;
    if (almC == 4) linSolveLM(); else linSolve();
''')
rep("    Ubase = U; real w = 1; int ls, lsOk = 0;\n", "    Ubase = U; lamBase = lamC; real w = 1; int ls, lsOk = 0;\n")
rep("        U = Ubase; dU = tv; dU *= w; U += dU; u1[] = U;\n", "        U = Ubase; dU = tv; dU *= w; U += dU; u1[] = U;\n        if (almC == 4) { lamC = dLam; lamC *= w; lamC += lamBase; }\n")
rep("        real rt = sqrt(rhs'*rhs);\n", "        real rt = sqrt(rhs'*rhs);\n        if (almC == 4) rt = sqrt(rt*rt + conRes());\n")
rep("      U = Ubase; u1[] = U; nLsFail++;\n", "      U = Ubase; u1[] = U; nLsFail++; if (almC == 4) lamC = lamBase;\n")
rep('''      du = w*tv.linfty/max(U.linfty, 1e-30);
      if (verbose) cout << "   it "''',
    '''      du = w*tv.linfty/max(U.linfty, 1e-30);
      if (kwcTan >= 2) {   // dual update: dw = (I - w g^T/n) dg/n + g/n - w, element step keeping |w| <= 1
        real[int] dgx = gsx - gox, dgy = gsy - goy, nr = gox .* gox, qq = goy .* goy; nr += qq; nr += dlt*dlt; nr = sqrt(nr);
        real[int] gd = gox .* dgx; qq = goy .* dgy; gd += qq; gd = gd ./ nr;          // g.dg/n
        real[int] ddx = dgx, ddy = dgy; qq = wdx .* gd; ddx -= qq; qq = wdy .* gd; ddy -= qq; ddx = ddx ./ nr; ddy = ddy ./ nr;
        qq = gox ./ nr; qq -= wdx; ddx += qq; qq = goy ./ nr; qq -= wdy; ddy += qq;
        for (int k = 0; k < nt; k++) {   // largest s in (0, 1] with |w + s dw| <= 1, times 0.99
          real A = ddx[k]^2 + ddy[k]^2, B = 2*(wdx[k]*ddx[k] + wdy[k]*ddy[k]), C = wdx[k]^2 + wdy[k]^2 - 1, sw = 1;
          if (A > 1e-300) { real disc = B*B - 4*A*C; if (disc > 0) { real sm = (-B + sqrt(disc))/(2*A); if (sm < 1) sw = 0.99*max(sm, 0.); } }
          wdx[k] += sw*ddx[k]; wdy[k] += sw*ddy[k];
        }
      }
      if (verbose) cout << "   it "''')

# ---- time step: multiplier set-up, cold-start iteration limit, retry, log ----
rep('''  tl = 0; ta0 = 0; ts = 0; tr = 0;
  it = newtonSolve(maxit, 0, "global it: "); du = duLast;
  bool okStep = (it < maxit) && (U.linfty < 1e30) && (U.linfty == U.linfty) && !nanHit;
''', '''  tl = 0; ta0 = 0; ts = 0; tr = 0;
  real wstep0 = wallNow();
  if (kwcTan >= 2) {
    if (!wInit) { kin(); real[int] nr = gsx .* gsx, qq = gsy .* gsy; nr += qq; nr += dlt*dlt; nr = sqrt(nr); wdx = gsx ./ nr; wdy = gsy ./ nr; wInit = true; }
    wsx = wdx; wsy = wdy;
  }
  if (almC) { lam0 = lamC; aMax = 1e30; }
  it = newtonSolve(maxit, 0, "global it: "); du = duLast;
  bool okStep = (it < maxit) && (U.linfty < 1e30) && (U.linfty == U.linfty) && !nanHit;
''')
rep('''      mark("retry: restore predictor + fresh tangent");
      localUpdate(); refreshTangent(); nretry++;
      it = newtonSolve(maxitRetry, 1000, "retry: "); du = duLast;
      okStep = (it < maxitRetry) && (U.linfty < 1e30) && (U.linfty == U.linfty) && !nanHit;
''', '''      mark("retry: restore predictor + fresh tangent");
      if (almC) { lamC = lam0; aMax = 1e30; }   // multiplier (a Newton unknown) and dual flux back to the step start
      if (kwcTan >= 2) { wdx = wsx; wdy = wsy; }
      localUpdate(); refreshTangent(); nretry++;
      it = newtonSolve(maxitRetry, 1000, "retry: "); du = duLast;
      okStep = (it < maxitRetry) && (U.linfty < 1e30) && (U.linfty == U.linfty) && !nanHit;
''')
rep('''  if (verbose) cout << "step " << n << " t " << t << " iters " << it << " |R| "''',
    '''  nlog << n << ",-1," << it << "," << (wallNow() - wstep0) << "," << aMax << "," << betaIt*180/pi << endl; nlog.flush; betaWorst = max(betaWorst, betaIt);
  if (verbose) cout << "step " << n << " t " << t << " iters " << it << " |R| "''')

# ---- -asm 1: direct assembly of the tangent matrix (plugin, fixed pattern) ----
rep('''  if (consReady && Jg.n == Wh.ndof) { matrix Kv = vKv(Wh, Wh); Aeff = Kconst + Kv; Aeff = Aeff + Jg; }
  else { matrix Kv = vKv(Wh, Wh); Aeff = Kconst + Kv; }
''', '''  if (asmFast) {
    asmChkNow = (asmCheck > 0) && (nAsmIt % max(asmCheck, 1) == 0); nAsmIt++;   // (FreeFEM evaluates both sides of &&)
    if (!asmReady) {   // once: pattern, entry order of the FreeFEM matrix, element -> entry map, constant part
      int[int] Ip(1), Jp(1); k2dAsmPattern(connA, Ip, Jp);
      { real[int] one(Ip.n); one = 1; Aeff = [Ip, Jp, one]; }
      //ASM-PREP
      int[int] I2(1), J2(1), Ik(1), Jk(1); real[int] C2(1), Ck(1);
      [I2, J2, C2] = Aeff; [Ik, Jk, Ck] = Kconst;
      real miss = k2dAsmSetup(connA, geoA, I2, J2, Ik, Jk, Ck);
      if (miss > 0 || I2.n != Ip.n) { cout << "ERROR: -asm 1: pattern mismatch (" << miss << " entries outside, " << I2.n << " vs " << Ip.n << ")" << endl; exit(1); }
      asmVals.resize(I2.n); asmX.resize(Wh.ndof); asmZ.resize(Wh.ndof); asmI.resize(I2.n); asmJ.resize(I2.n); asmI = I2; asmJ = J2;
      { real[int] dg(Wh.ndof); dg = Kconst.diag;   // probe vectors: zero on the prescribed dofs (their diagonal carries tgv)
        for (int i = 0; i < Wh.ndof; i++) { real fr = (abs(dg[i]) < 1e20) ? 1. : 0.; asmX[i] = sin(1. + 0.7*i)*fr; asmZ[i] = cos(0.3*i)*fr; } }
      asmReady = true;
    }
    asmHK(:, 0) = H11[]; asmHK(:, 1) = H12[]; asmHK(:, 2) = H22[]; asmHK(:, 3) = Hd[];
    k2dAsmA(asmHK, ((consReady && asmJok) ? 1 : 0), asmVals);
    if (asmSerialDirect && (linMode == 0 || nLinFall != asmFallSeen)) { Aeff = [asmI, asmJ, asmVals]; asmFallSeen = nLinFall; }
    else Aeff.coef = asmVals;
    { // the entry order of the matrix must be the one the map was built for: compare z'Ax of the matrix with the plugin's
      // (scale: the sum of |z_i A_ij x_j|, since z'Ax itself can pass through zero)
      real[int] yp = Aeff*asmX; yp = yp .* asmZ; real s1 = yp.sum, s2 = k2dAsmProbe(asmVals, asmX, asmZ), sa = k2dAsmProbeAbs(asmVals, asmX, asmZ);
      asmDev = max(asmDev, abs(s1 - s2)/max(sa, 1e-290));
      if (!(abs(s1 - s2) <= 1e-9*sa + 1e-290)) {
        cerr << "ERROR: -asm 1: the assembled matrix does not match the plugin's (z'Ax " << s1 << " vs " << s2 << ", scale " << sa << "): entry order changed, or not-a-number in the tangent" << endl; exit(1); }
    }
    if (asmChkNow && jgFresh) {   // ASMCHECK: the same matrix by the forms, from the same element blocks (rank-local in MPI)
      real tck0 = clock();
      matrix Kv = vKv(Wh, Wh); matrix As = Kconst + Kv; if (consReady && asmJok) As = As + Jg;
      int[int] ci(1), cj(1); real[int] cc(1); [ci, cj, cc] = As;
      real dev = gmaxS(k2dAsmCmp(ci, cj, cc, asmVals)); nAsmChk++; asmChkMax = max(asmChkMax, dev);
      if (verbose) cout << "   ASMCHECK " << nAsmChk << " (Newton iteration " << nAsmIt << " of the run): max |direct - forms| / max |A| = " << dev << " (" << (clock() - tck0) << " s)" << endl;
      if (!(dev <= 1e-12)) {
        if (isRoot()) cout << "ERROR: ASMCHECK at Newton iteration " << nAsmIt << " of the run: the direct assembly differs from the forms by " << dev << " of the largest entry (limit 1e-12); stopping. Restart with FASTASM=0 to continue with the form assembly" << endl;
        exit(1);
      }
    }
    jgFresh = false;
  }
  else if (consReady && Jg.n == Wh.ndof) { matrix Kv = vKv(Wh, Wh); Aeff = Kconst + Kv; Aeff = Aeff + Jg; }
  else { matrix Kv = vKv(Wh, Wh); Aeff = Kconst + Kv; }
''')
rep('''  real[int] tv(3*nt);
  packPd(HQ, 0, tv);''', '''  if (asmFast) { k2dAsmJ(HQ, CQ); asmJok = true; if (!(asmCheck > 0 && nAsmIt % max(asmCheck, 1) == 0)) return 0; jgFresh = true; }   // ASMCHECK iteration: the form blocks too
  real[int] tv(3*nt);
  packPd(HQ, 0, tv);''')

rep('cout << "done. wall time "', 'if (asmFast) cout << "fast assembly: largest matrix check deviation " << asmDev << " (limit 1e-9; rank 0)" << endl;\nif (betaMax > 0) cout << "branch guard (-betamax " << betaMax*180/pi << " deg): largest |beta| at accepted steps " << betaWorst*180/pi << " deg; rejected iterates " << nBetaIt << ", rejected steps " << nBetaRej << endl;\nif (locTol > 0 || betaMax > 0) cout << "final-state checks: rejected for the global residual, the constraint or admissibility " << nFinRej << endl;\nif (locTol > 0) cout << "local acceptance test (-localtol " << locTol << "): largest residual " << locWorst << ", rejections " << nLocRej << "; slip law " << (slipExact ? "exact (no cap)" : "capped at x = 300") << endl;\nif (asmCheck > 0) cout << "ASMCHECK: " << nAsmChk << " comparisons with the form assembly, largest deviation " << asmChkMax << " of the largest entry (limit 1e-12)" << endl;\ncout << "done. wall time "')

# ---- whole-step timing (-verbose 1): a second line per step with everything outside the Newton timers ----
rep("    if (!lsOk) {\n", "    tLs += clock() - tLs0;\n    if (!lsOk) {\n")
rep("    Ubase = U; lamBase = lamC; real w = 1; int ls, lsOk = 0;\n", "    Ubase = U; lamBase = lamC; real w = 1; int ls, lsOk = 0; real tLs0 = clock();\n")
rep("for (int n = nstart + 1; n <= ntot; n++) {\n  bool relax = (n <= nrelax);\n", "for (int n = nstart + 1; n <= ntot; n++) {\n  bool relax = (n <= nrelax);\n  tStep0 = clock(); tLs = 0;\n")
rep("  tl = 0; ta0 = 0; ts = 0; tr = 0;", "  tl = 0; ta0 = 0; ts = 0; tr = 0; tNewt0 = clock();")
rep('''  if (verbose) cout << "step " << n << " t " << t << " iters "''', '''  tNewt1 = clock();
  if (verbose) cout << "step " << n << " t " << t << " iters "''')
rep('''  if (clock() - clk0 > walllimit) { saveCk(n);''', '''  if (verbose) cout << "step " << n << " time (rank 0, CPU s): total " << (clock() - tStep0) << " = before Newton " << (tNewt0 - tStep0) << " + Newton " << (tNewt1 - tNewt0)
    << " (of it line search " << tLs << ") + after Newton (phase field, temperature, densities, output, checkpoints) " << (clock() - tNewt1) << endl;
  if (clock() - clk0 > walllimit) { saveCk(n);''')

# ---- exact slip rate law (-slipexact 1): no cap on the stress ratio, derivative of the same expression ----
rep('''func int flow(real[int] & t, real f, real[int] & gd, real[int] & dgd, real[int] & eo, int rel) {
  effDrag(t, rel); eo = etaE;
''', '''// v2 (-slipexact 1): 1/|gd| = 1/(g0 x^m) + eE/|t| evaluated without overflow for every x = |t|/(f tau_c).  With L = ln(g0 eE/(f tau_c)) + (m-1) ln x
// (the log of the ratio drag term / activation term) and sigma = 1/(1 + e^-L):  |gd| = sigma |t|/eE,
// d|gd|/d|t| = (sigma/eE) [m (1 - sigma) + sigma] - (sigma/eE) |gd| d eE/d|t|.  d ln|gd| / d ln x lies in [1, m] at fixed eE and falls towards 0 with the relativistic drag near the limiting rate.
int slipExact = atoi(arg("-slipexact", "0"));
// -betamax b (radians): lattice-rotation branch guard. beta = polar angle of Fe = F Fp^-1 minus theta = atan2(Ub21 - Ub12, Ub11 + Ub22)
// (Ubar = R(beta) S, S > 0); the Cosserat energy is pi-periodic in beta (a = tr(S) sin(beta)/2 also vanishes at beta = pi). A trial
// iterate with |beta| > 2b anywhere is inadmissible (line-search backtracking, like the det F bound); an accepted step needs |beta| <= b.
real betaMax = atof(arg("-betamax", "0")), betaIt = 0, betaWorst = 0, betaQ = 0; int nBetaIt = 0, nBetaRej = 0;
func int flowExact(real[int] & t, real f, real[int] & gd, real[int] & dgd, int rel) {
  real[int] at = abs(t), LL(nt), sg(nt), rr(nt), qq(nt);
  xa = at ./ tauc; xa *= 1./f; tmp = xa; tmp += 1e-300; lx = log(tmp);
  LL = etaE ./ tauc; LL *= g0/f; LL = log(LL); tmp = lx; tmp *= (mexp - 1); LL += tmp;
  sg = -LL; sg = exp(sg); sg += 1; sg = one ./ sg;              // sigma (exp overflow -> sigma = 0: activation-limited)
  rr = exp(LL); rr += 1; rr = one ./ rr;                        // 1 - sigma without cancellation
  qq = sg ./ etaE;                                              // |gd|/|t|
  gd = at .* qq;
  tmp = rr; tmp *= mexp; tmp += sg; dgd = qq .* tmp;
  if (rel && vcrel > 0) { tmp = qq .* gd; tmp = tmp .* detaE; dgd -= tmp; }
  return 0;
}
func int flow(real[int] & t, real f, real[int] & gd, real[int] & dgd, real[int] & eo, int rel) {
  effDrag(t, rel); eo = etaE;
  if (slipExact) { flowExact(t, f, gd, dgd, rel); return 0; }
''')

# ---- branch measure at the three elastic Gauss points (in gaussInt) and at the element means ----
rep("func int gaussInt(int wantF) {\n", """func real betaNow() {   // max over the elements of |beta| = |atan2(Ub21 - Ub12, Ub11 + Ub22)| for the current Ubar arrays (rank-local)
  real[int] cb = ub11 + ub22, sb = ub21 - ub12;
  if (cb.min > 0) { sb = abs(sb); sb = sb ./ cb; return atan(sb.max); }   // all |beta| < pi/2: vectorised
  real bm = 0; for (int k = 0; k < cb.n; k++) bm = max(bm, abs(atan2(sb[k], cb[k]))); return bm;
}
func int gaussInt(int wantF) {
""")
rep("    if (almC) { if (q == 0) aMean = 0; tmp = aa; tmp *= gqw(q); aMean += tmp; }\n",
    "    if (almC) { if (q == 0) aMean = 0; tmp = aa; tmp *= gqw(q); aMean += tmp; }\n    if (betaMax > 0) { if (q == 0) betaQ = 0; betaQ = max(betaQ, betaNow()); }   // branch measure at the Gauss point\n")

# ---- branch guard (-betamax): every residual evaluation measures |beta| at the element-mean state ----
rep('''  gaussInt(1);
  rhs = rin;
''', '''  gaussInt(1);
  if (betaMax > 0) {   // the three elastic Gauss points (betaQ, from gaussInt) and the element means (gaussInt ends there)
    betaIt = gmaxS(max(betaQ, betaNow()));
    if (betaIt > 2*betaMax) { if (!jBad) nBetaIt++; jBad = 1; jbLast = -1e30; if (verbose) cout << "   inadmissible iterate: lattice-rotation mismatch |beta| = " << betaIt*180/pi << " deg > " << 2*betaMax*180/pi << " deg" << endl; }
  }
  rhs = rin;
''')

# ---- the det F counter also counts the iterates the branch guard rejects: report them apart ----
rep('''cout << "safeguard counters: det F bound rejections " << nJBad << ", ''', '''cout << "safeguard counters: det F bound rejections " << max(nJBad - nBetaIt, 0) << ", branch-guard rejections " << nBetaIt << ", ''')

# ---- local acceptance test (-localtol r): the local constitutive equations must hold at the accepted state ----
rep("func int newtonSolve(int mx, int itOff, string tg) {\n", '''real locTol = atof(arg("-localtol", "0")); int nLocRej = 0; real locWorst = 0;   // -localtol: limit on the local residuals (increment units)
func real localAcceptRes() { real[int] rn(nt); localRes(rn); return gmaxS(rn.max); }   // largest |r_1| + |r_2| + |r_3| + |r_*| over the elements
// final state of a converged step (with -localtol or -betamax): ONE last local update, then the global residual (and the reaction
// on the prescribed dofs), the constraint, the local residuals and the branch measure, all of that same state, which is the state
// committed (the final local update of the time loop is skipped: stateFinal). Returns 1 if the state is accepted.
int stateFinal = 0, nFinRej = 0;
func int acceptFinal(int n, string pf) {
  localUpdate(); residual(); stateFinal = 1;
  real mf = sqrt(rhs'*rhs); if (almC == 4) mf = sqrt(mf*mf + conRes());
  real am = 0; if (almC) { real[int] qa = abs(aMean); am = gmaxS(qa.max); aMax = am; }
  int ok = 1;
  if (jBad) { ok = 0; nFinRej++; cout << pf << "final state inadmissible (det F bound, |q| > 50 or |beta| > 2 betamax)" << endl; }
  if (ok && !(mf <= max(resLast*(1 + 1e-3), tolR*max(res0Last, rFloor)))) { ok = 0; nFinRej++; cout << pf << "global residual of the final state " << mf << " above the accepted " << resLast << endl; }
  if (ok && almC && !(am <= almFin*almTol)) { ok = 0; nFinRej++; cout << pf << "constraint of the final state max |abar| " << am << " above " << almFin*almTol << " (-almfinal x -almtol)" << endl; }
  if (ok && locTol > 0) {
    real lr = localAcceptRes(); locWorst = max(locWorst, lr);
    if (!(lr <= locTol)) { ok = 0; nLocRej++; cout << pf << "local constitutive residual " << lr << " above -localtol " << locTol << " at the final state" << endl; }
  }
  if (ok && betaMax > 0 && !(betaIt <= betaMax)) { ok = 0; nBetaRej++; cout << pf << "lattice-rotation mismatch |beta| = " << betaIt*180/pi << " deg above -betamax " << betaMax*180/pi << " deg at the final state (another branch of the Cosserat energy)" << endl; }
  return ok;
}
func int newtonSolve(int mx, int itOff, string tg) {
''')
rep('''  bool okStep = (it < maxit) && (U.linfty < 1e30) && (U.linfty == U.linfty) && !nanHit;
''', '''  bool okStep = (it < maxit) && (U.linfty < 1e30) && (U.linfty == U.linfty) && !nanHit;
  stateFinal = 0;
  if (okStep && (locTol > 0 || betaMax > 0)) okStep = acceptFinal(n, "step " + n + ": ");   // rejected: the retry path below
''')
rep('''      okStep = (it < maxitRetry) && (U.linfty < 1e30) && (U.linfty == U.linfty) && !nanHit;
''', '''      okStep = (it < maxitRetry) && (U.linfty < 1e30) && (U.linfty == U.linfty) && !nanHit;
      stateFinal = 0;
      if (okStep && (locTol > 0 || betaMax > 0)) okStep = acceptFinal(n, "   retry: ");
''')

# ---- the state accepted by acceptFinal is the one committed: no further local update ----
rep('''  mark("converged: final local update"); localUpdate();
''', '''  mark("converged: final local update"); if (!stateFinal) localUpdate(); stateFinal = 0;   // v2: acceptFinal made (and checked) the final state
''')

# ---- spatial power diagnostic: centroid driving forces (the local law) against their Gauss-point mean (conjugate to the energy) ----
rep("func int newtonSolve(int mx, int itOff, string tg) {\n", """// The plastic state is P0 and evolves with the driving forces at the element mean (centroid) state; the elastic energy is sampled at
// the three Gauss points, whose conjugate driving force for a uniform element increment is the Gauss mean of tau. powerSplit()
// evaluates both at the accepted state: P_c = sum_e A_e tau_c . dq, P_q = sum_e A_e <tau>_q . dq (dq: the three slips and dw*).
real powCcum = 0, powQcum = 0, powAbsCum = 0, powLocCum = 0, powDevMax = 0, powLocMax = 0, tauDevMax = 0; int powN0 = -1, powN1 = -1;
// powerSplit() evaluates at the accepted mechanical state and buffers; powerPublish() writes the row and adds to the totals only once
// the whole coupled step has succeeded (phase field, densities, heat, commit without NaN). Per step, rank 0: <out>_power.csv;
// run.sh drops rows after the checkpoint on a resume, so the file holds the active history only
int pbN = -1, pbAcc = 0; real pbT = 0, pbC = 0, pbQ = 0, pbA = 0, pbL = 0, pbD = 0;
string plogFile = isRoot() ? (out + "_power.csv") : "/dev/null"; bool plHdr = (exec("test -s " + plogFile) != 0);
ofstream plog(plogFile, append); if (plHdr) plog << "step,t_ns,P_c,P_q,sum_abs_P_c,sum_abs_dP,max_rel_dtau" << endl;
func int powerSplit(int nstep, real tnow, bool acc) {   // acc: add to the totals (loaded steps only)
  real[int] q1(nt), q2(nt), q3(nt), qs(nt), w(nt); q1 = 0; q2 = 0; q3 = 0; qs = 0;
  real[int] c1 = ta, c2 = tb, c3 = tcc, cs = tstar;   // centroid values of the accepted state (taus() of the last local update)
  trialFp();
  for (int q = 0; q < nQ; q++) {
    e11 = GD1[q]*U; hy1 = GD2[q]*U; hx2 = GD3[q]*U; e22 = GD4[q]*U; thc = GD5[q]*U;
    finiteStress();
    w = Mn11 .* A11a; tmp = Mn12 .* A12a; w += tmp; tmp = Mn21 .* A21a; w += tmp; tmp = Mn22 .* A22a; w += tmp; w *= gqw(q); q1 += w;
    w = Mn11 .* A11b; tmp = Mn12 .* A12b; w += tmp; tmp = Mn21 .* A21b; w += tmp; tmp = Mn22 .* A22b; w += tmp; w *= gqw(q); q2 += w;
    w = Mn11 .* A11c; tmp = Mn12 .* A12c; w += tmp; tmp = Mn21 .* A21c; w += tmp; tmp = Mn22 .* A22c; w += tmp; w *= gqw(q); q3 += w;
    w = Mn21 - Mn12; w *= gqw(q); qs += w;
  }
  kin(); taus();   // back to the centroid state (and its tau) used by the commit
  real[int] pc = c1 .* dpa, pq = q1 .* dpa, pa(nt);
  tmp = c2 .* dpb; pc += tmp; tmp = c3 .* dpc; pc += tmp; tmp = cs .* dwst; pc += tmp;
  tmp = q2 .* dpb; pq += tmp; tmp = q3 .* dpc; pq += tmp; tmp = qs .* dwst; pq += tmp;
  pa = abs(pc); pa = pa .* areas; real sa = gsumEs(pa);
  tmp = pc .* areas; real sc = gsumEs(tmp); tmp = pq .* areas; real sq = gsumEs(tmp);
  w = pq - pc; w = abs(w); w = w .* areas; real sl = gsumEs(w);   // non-cancelling: element by element
  pbN = nstep; pbT = tnow; pbAcc = acc; pbC = sc; pbQ = sq; pbA = sa; pbL = sl;
  real tm = max(max(c1.linfty, c2.linfty), max(c3.linfty, cs.linfty));
  w = q1 - c1; real dm = w.linfty; w = q2 - c2; dm = max(dm, w.linfty); w = q3 - c3; dm = max(dm, w.linfty); w = qs - cs; dm = max(dm, w.linfty);
  pbD = gmaxS(dm)/max(gmaxS(tm), 1e-300);
  return 0;
}
func int powerPublish(int nstep) {   // the step succeeded: publish its buffered values
  if (pbN != nstep) return 0;
  plog.precision(10); plog << pbN << "," << pbT << "," << pbC << "," << pbQ << "," << pbA << "," << pbL << "," << pbD << endl; plog.flush;
  if (pbAcc) {
    powCcum += pbC; powQcum += pbQ; powAbsCum += pbA; powLocCum += pbL; if (powN0 < 0) powN0 = pbN; powN1 = pbN;
    if (pbA > 0) { powDevMax = max(powDevMax, abs(pbQ - pbC)/pbA); powLocMax = max(powLocMax, pbL/pbA); }
    tauDevMax = max(tauDevMax, pbD);
  }
  pbN = -1; return 0;
}
func int newtonSolve(int mx, int itOff, string tg) {
""")
rep('  mark("commit: Newmark update, external work, dissipation");\n', '  powerSplit(n, t, !relax);\n  mark("commit: Newmark update, external work, dissipation");\n')
rep('cout << "done. wall time "', 'cout << "spatial power split, loaded steps " << powN0 << ".." << powN1 << " of this launch (every step of the run: <out>_power.csv): plastic work with centroid driving forces " << powCcum << ", with their Gauss-point mean " << powQcum << "; net difference " << abs(powQcum - powCcum)/max(powAbsCum, 1e-300) << " and element-by-element sum |P_q - P_c| " << powLocCum/max(powAbsCum, 1e-300) << " of sum |P_c| (largest per step " << powDevMax << ", " << powLocMax << "); largest |<tau>_q - tau_c| / max |tau_c| " << tauDevMax << endl;\ncout << "done. wall time "')

# ---- publish the power-split values of a step once it has succeeded (after the NaN check of the commit) ----
rep('  mark("output: checkpoints");\n', '  powerPublish(n);\n  mark("output: checkpoints");\n')

# ---- failed step: the checkpoint of step n-1 gets the multiplier and dual flux of step n-1 (not those of the failed attempt) ----
rep('''    if (!okStep) {   // stop cleanly at the end of step n-1 (nothing of step n is committed yet except U)
      U = Un; u1[] = U; dpa = 0; dpb = 0; dpc = 0; dwst = 0;
''', '''    if (!okStep) {   // stop cleanly at the end of step n-1 (nothing of step n is committed yet except U)
      U = Un; u1[] = U; dpa = 0; dpb = 0; dpc = 0; dwst = 0;
      if (almC) lamC = lam0;                       // v2 state at the end of step n-1 for the checkpoint below
      if (kwcTan >= 2) { wdx = wsx; wdy = wsy; }
''')
# ---- checkpoint publication requires the .mult side file when it carries state ----
rep('''    chk = chk + " && test -s " + gen + ".gnd";
''', '''    chk = chk + " && test -s " + gen + ".gnd";
    if (almC || kwcTan >= 2) chk = chk + " && test -s " + gen + ".mult";
''')
# ---- VTU psi_elastic: the physical elastic energy density (without the multiplier term lambda a) ----
rep('''f29[] = psiF;''', '''f29[] = psiF; if (almC) { real[int] la = lamC .* aa; f29[] -= la; } /* v2: psi_elastic without lambda a (physical energy) */''')

# ---- checkpoints: multiplier lambda and dual flux w in the side file <checkpoint>.mult ----
rep('''  { ofstream f(fn + ".gnd"); f.precision(17); f << n << endl << GR1 << GR2; }   // Curl-Fp GND reference (side file)
''', '''  { ofstream f(fn + ".gnd"); f.precision(17); f << n << endl << GR1 << GR2; }   // Curl-Fp GND reference (side file)
  if (almC || kwcTan >= 2) { if (!wInit) { wdx = 0; wdy = 0; }
    { ofstream f(fn + ".mult"); f.precision(17); f << n << endl << (wInit ? 1 : 0) << endl << lamC << wdx << wdy; } }   //MULT-SAVE (Cosserat multiplier, KWC dual flux)
''')
rep("for x in '' .rhoM .gnd .plate; do", "for x in '' .rhoM .gnd .plate .mult; do")
rep('''  u1[] = U;
  { // mobile density: side file <restart>.rhoM''', '''  u1[] = U;
  if (almC || kwcTan >= 2) {   // multiplier and dual flux (absent: rebuilt in the first step)
    int nmu = -1, wi = 0;
    if (exec("test -s " + restart + ".mult") == 0) { ifstream fmu(restart + ".mult"); fmu >> nmu; if (nmu == nn) { fmu >> wi;
      fmu >> lamC >> wdx >> wdy;   //MULT-READ
      multRead = true; } }
    if (multRead) { wInit = (wi == 1); cout << "multiplier / dual flux read from " << restart << ".mult" << endl; }
    else cout << "no " << restart << ".mult for step " << nn << ": multiplier / dual flux start cold" << endl;
  }
  { // mobile density: side file <restart>.rhoM''')

# ---- metadata ----
rep('''  f << " \\"tangent\\": {''',
    '''  f << " \\"v2\\": {\\"kwctan\\": " << kwcTan << ", \\"almc\\": " << almC << ", \\"asm\\": " << asmFast << ", \\"slipexact\\": " << slipExact << ", \\"localtol\\": " << locTol << ", \\"betamax_rad\\": " << betaMax << "},\\n";
  f << " \\"tangent\\": {''')
# ---- patch_retrytan_r7.py: second retry with the exact KWC Hessian ----
def _rt_rep(a, b):
    global s
    n = s.count(a); assert n == 1, ('retrytan patch target not found once', a[:80], n); s = s.replace(a, b)
_rt_rep('int maxitRetry = atoi(arg("-maxitretry", "80"));',
        'int maxitRetry = atoi(arg("-maxitretry", "80"));\nint retryTan = atoi(arg("-retrytan", "-1")), maxitRetry2 = atoi(arg("-maxitretry2", "80")), nRetryTan = 0, nRetryTanOk = 0;   // second retry with another KWC tangent (patch_retrytan_r7.py)')
_rt_a = '"   retry with fresh tangent: " << (okStep ? "converged" : "FAILED") << " in " << it << " its" << endl;'
assert s.count(_rt_a) == 1, 'retrytan: retry line not found once'
_i = s.index(_rt_a) + len(_rt_a); _i = s.index('\n', _i - 1) + 1 if s[_i:_i+1] != '\n' else _i + 1
s = s[:_i] + """      if (!okStep && !nanHit && kwcTan >= 2 && retryTan >= 0 && retryTan != kwcTan) {   // second retry: same problem, exact KWC Hessian (patch_retrytan_r7.py)
        U = Upred; u1[] = U; dpa = 0; dpb = 0; dpc = 0; dwst = 0;
        if (almC) { lamC = lam0; aMax = 1e30; }
        wdx = wsx; wdy = wsy;
        int kt0 = kwcTan; kwcTan = retryTan; nRetryTan++;
        mark("retry 2: restore predictor, KWC tangent -kwctan " + retryTan);
        localUpdate(); refreshTangent(); nretry++;
        it = newtonSolve(maxitRetry2, 2000, "retry 2: "); du = duLast;
        real uL2 = gmaxS(U.linfty); okStep = (it < maxitRetry2) && (uL2 < 1e30) && (uL2 == uL2) && !nanHit;
        stateFinal = 0;
        if (okStep && (locTol > 0 || betaMax > 0)) okStep = acceptFinal(n, "   retry 2: ");
        kwcTan = kt0;
        if (okStep) { nRetryTanOk++; kin(); real[int] nr = gsx .* gsx, qq = gsy .* gsy; nr += qq; nr += dlt*dlt; nr = sqrt(nr); wdx = gsx ./ nr; wdy = gsy ./ nr; }   // dual flux of the converged state
        cout << "   retry with the exact KWC Hessian (-kwctan " << retryTan << "): " << (okStep ? "converged" : "FAILED") << " in " << it << " its" << endl;
      }
""" + s[_i:]
_rt_d = 'cout << "done. wall time '
assert s.count(_rt_d) == 1, 'retrytan: end-of-run line not found once'
_j = s.rfind('\n', 0, s.index(_rt_d)) + 1
s = s[:_j] + 'cout << "second retries with the exact KWC Hessian (-retrytan " << retryTan << "): " << nRetryTan << ", converged " << nRetryTanOk << endl;   // patch_retrytan_r7.py\n' + s[_j:]
# ---- patch_localcold_r7.py: local update from zero increments at every call ----
def _lc_rep(a, b):
    global s
    n = s.count(a); assert n == 1, ('localcold patch target not found once', a[:80], n); s = s.replace(a, b)
_lc_rep('real locres = 0; int nlocit = 0;',
        'real locres = 0; int nlocit = 0;\nint locCold = atoi(arg("-localcold", "1")); real locResMax = 0;   // patch_localcold_r7.py: local update from zero increments; largest local residual left')
_lc_rep('func int localUpdate() {\n  kin();\n',
        'func int localUpdate() {\n  kin();\n  if (locCold) { dpa = 0; dpb = 0; dpc = 0; dwst = 0; }   // patch_localcold_r7.py: the local state is a function of (U, lambda) only\n')
_lc_rep('  polishLocal();\n  taus();\n  return 0;\n}',
        '  polishLocal();\n  if (locCold) { real[int] rq(nt); localRes(rq); locResMax = max(locResMax, rq.max); }   // also makes the slip rates consistent with the final increments\n  taus();\n  return 0;\n}')
_lc_d = 'cout << "done. wall time '
assert s.count(_lc_d) == 1, 'localcold: end-of-run line not found once'
_lc_j = s.rfind('\n', 0, s.index(_lc_d)) + 1
s = s[:_lc_j] + '{ real lrm = gmaxS(locResMax); cout << "local update (-localcold " << locCold << "): largest local residual left by any call " << lrm << " (increment units)" << endl; }   // patch_localcold_r7.py\n' + s[_lc_j:]
# ---- patch_localsolve_r7.py: local polish to convergence per element; warm start with a cold restart of unconverged elements ----
def _ls_rep(a, b):
    global s
    n = s.count(a); assert n == 1, ('localsolve patch target not found once', a[:80], n); s = s.replace(a, b)
_ls_rep('arg("-localcold", "1")', 'arg("-localcold", "0")')
_ls_rep('real locres = 0; int nlocit = 0;', 'real locres = 0; int nlocit = 0;\nint locMaxIt = atoi(arg("-localmaxit", "200")), nPolishIt = 0, nLocRestart = 0;   // patch_localsolve_r7.py')
_ls_rep('  if (locCold) { dpa = 0; dpb = 0; dpc = 0; dwst = 0; }   // patch_localcold_r7.py: the local state is a function of (U, lambda) only\n', '')
_ls_rep('  if (locCold) { real[int] rq(nt); localRes(rq); locResMax = max(locResMax, rq.max); }   // also makes the slip rates consistent with the final increments\n', '')
_ls_rep('func int localUpdate() {\n  kin();\n', 'func int localSolve() {   // one local solve from the current increments (patch_localsolve_r7.py)\n  kin();\n')
_ls_end = '  polishLocal();\n  taus();\n  return 0;\n}\n'
_ls_i = s.index('func int localSolve() {'); _ls_j = s.index(_ls_end, _ls_i) + len(_ls_end)
s = s[:_ls_j] + """func int localUpdate() {   // patch_localsolve_r7.py
  if (locCold) { dpa = 0; dpb = 0; dpc = 0; dwst = 0; }
  localSolve();
  real[int] rq(nt); localRes(rq);
  real[int] kp(nt); for [i, v : kp] v = (rq[i] > 1e-12);   // 1 where not converged
  real nbad = gmaxS(kp.max);   // the same on every rank, so that all ranks take the same branch
  if (nbad > 0.5) {   // elements not converged from the warm start: solve those once more from zero increments
    kp = one - kp;
    dpa = dpa .* kp; dpb = dpb .* kp; dpc = dpc .* kp; dwst = dwst .* kp; nLocRestart++;
    localSolve(); localRes(rq);
  }
  locResMax = max(locResMax, rq.max);
  taus();
  return 0;
}
""" + s[_ls_j:]
_ls_pi = s.index('func int polishLocal() {'); _ls_pj = s.index('\n}\n', _ls_pi) + 3
s = s[:_ls_pi] + r"""func int polishLocal() {   // coupled slip + spin Newton with the exact local Jacobian, per element; runs until every element has converged
  // (sum |r| < 1e-14) or no unconverged element can still lower its residual. Per-element backtracking 1, 1/2, ..., 2^-11;
  // converged elements are left untouched. (patch_localsolve_r7.py; before: at most -polish iterations, stopped as soon as the
  // largest residual over all elements fell by less than half, which left slow elements unconverged)
  real[int] rn0(nt), rn1(nt), qa(nt), qb(nt), qc(nt), qw(nt), ok(nt), pend(nt), acc(nt), act(nt), da(nt), db(nt), dc(nt), dw(nt);
  int cap = max(nPolish, locMaxIt);
  for (int ip = 0; ip < cap; ip++) {
    localRes(rn0); lastPolishRes = rn0.max; nPolishIt = ip;
    if (lastPolishRes < 1e-14) break;
    act = rn0; act -= 1e-14; tmp = abs(act); act += tmp; tmp *= 2; tmp += 1e-300; act = act ./ tmp;   // 1 where sum |r| > 1e-14
    elemJac(); buildJqq();
    tmp = -ra; BB(:, 0) = tmp; tmp = -rb; BB(:, 1) = tmp; tmp = -rc; BB(:, 2) = tmp; tmp = -rwl; BB(:, 3) = tmp;
    solve4(1);
    da = BB(:, 0); db = BB(:, 1); dc = BB(:, 2); dw = BB(:, 3);
    qa = dpa; qb = dpb; qc = dpc; qw = dwst;
    pend = act; acc = one - act;   // converged elements count as accepted with their current value
    real w = 1; real[int] prog(nt); prog = 0;
    for (int ls = 0; ls < 12; ls++) {
      tmp = da; tmp *= w; tmp += qa; tmp = tmp .* pend; tmp2 = dpa .* acc; dpa = tmp + tmp2;
      tmp = db; tmp *= w; tmp += qb; tmp = tmp .* pend; tmp2 = dpb .* acc; dpb = tmp + tmp2;
      tmp = dc; tmp *= w; tmp += qc; tmp = tmp .* pend; tmp2 = dpc .* acc; dpc = tmp + tmp2;
      tmp = dw; tmp *= w; tmp += qw; tmp = tmp .* pend; tmp2 = dwst .* acc; dwst = tmp + tmp2;
      if (spinCap > 0) { real[int] smx(nt); smx = spinCap; clampv(dwst, smx); }   // principal branch also in the polish
      localRes(rn1);
      tmp = rn0 - rn1; tmp2 = abs(tmp); tmp += tmp2; tmp2 *= 2; tmp2 += 1e-300; ok = tmp ./ tmp2;   // 1 where the residual decreased
      ok = ok .* pend; acc += ok; prog += ok; pend -= ok;
      if (pend.max < 0.5) break;
      w *= 0.5;
    }
    // still pending: back to the old increments
    tmp = dpa .* acc; tmp2 = qa .* pend; dpa = tmp + tmp2; tmp = dpb .* acc; tmp2 = qb .* pend; dpb = tmp + tmp2;
    tmp = dpc .* acc; tmp2 = qc .* pend; dpc = tmp + tmp2; tmp = dwst .* acc; tmp2 = qw .* pend; dwst = tmp + tmp2;
    if (gmaxS(prog.max) < 0.5) break;   // (all ranks together) no unconverged element could lower its residual: round-off floor (or a stuck element)
  }
  return 0;
}
""" + s[_ls_pj:]
_ls_rep('largest local residual left by any call " << lrm << " (increment units)" << endl;',
        'largest local residual left by any call " << lrm << " (increment units), local updates restarted from zero increments " << nLocRestart << endl;')
# ---- patch_lsarmijo_r7.py: Armijo constant of the global line search ----
def _la_rep(a, b):
    global s
    n = s.count(a); assert n == 1, ('lsarmijo patch target not found once', a[:80], n); s = s.replace(a, b)
_la_rep('int lsMax = atoi(arg("-lsmax", "6"));', 'int lsMax = atoi(arg("-lsmax", "6"));\nreal lsArm = atof(arg("-lsarmijo", "0.25"));   // Armijo constant of the line search on the merit (patch_lsarmijo_r7.py; was 1e-4)')
_la_rep('rt <= (1 - 1e-4*w)*rn', 'rt <= (1 - lsArm*w)*rn')
# ---- patch_gndcg_r7.py: Cermelli-Gurtin GND (planar vector c/Jp), no reference subtraction ----
def _cg_rep(a, b):
    global s
    n = s.count(a); assert n == 1, ('gndcg patch target not found once', a[:80], n); s = s.replace(a, b)
_cg_rep('int gndCG = 0;   // curl convention recorded in the checkpoint manifest (0: G = Fp c)',
        'int gndCG = 1;   // 1: Cermelli-Gurtin g = c/Jp, rho_G = |g|/b, no reference (patch_gndcg_r7.py); recorded in the checkpoint manifest\n'
        'cout << "GND measure: Cermelli-Gurtin g = c/Jp (planar curl of Fp), rho_G = |g|/b, no reference subtraction" << endl;')
_cg_rep('  g1 = fp11 .* c1e[]; tmp = fp12 .* c2e[]; g1 += tmp;   // G = Fp c\n  g2 = fp21 .* c1e[]; tmp = fp22 .* c2e[]; g2 += tmp;',
        '  {   // Cermelli-Gurtin, strict plane strain: G_CG = e3 (x) (c/Jp); stored vector g = c/Jp (Jp = det Fp = 1 for the exponential update)\n'
        '    real[int] jpd = fp11 .* fp22; tmp = fp12 .* fp21; jpd -= tmp; for [i, v : jpd] v = max(v, 1e-6);\n'
        '    g1 = c1e[] ./ jpd; g2 = c2e[] ./ jpd;\n'
        '  }')
# restart: the stored reference of an old checkpoint is ignored
_cg_rep('if (ng == nn && max(manGndcg, 0) != gndCG) {   // the stored reference was computed with the other curl convention: never subtract it\n'
        '      cout << "RESTART ERROR: the stored GND reference uses -gndcg " << max(manGndcg, 0) << " but this run has -gndcg " << gndCG << ": restart with the same -gndcg" << endl; exit(4); }',
        'GR1 = 0; GR2 = 0;   // the stored reference of an old checkpoint is read and discarded\n    if (ng == nn && max(manGndcg, 0) != gndCG) cout << "NOTE: the checkpoint was written with the old GND measure (G = Fp c with a relaxation reference); rho_G is now |c/Jp|/b from the current Fp, the stored reference is ignored" << endl;')
_cg_rep('    if (ng == nn) cout << "GND reference read from " << restart << ".gnd" << endl;\n', '')
_cg_rep('    else if (nn != nrelax) cout << "WARNING: no " << restart << ".gnd for step " << nn << ": GND reference set to 0" << endl;\n', '')
# rho_G from the current Fp at all times, including the state at the end of the relaxation
_cg_rep('curlFpG(g1, g2); g1 -= GR1; g2 -= GR2; tmp = g1 .* g1;', 'curlFpG(g1, g2); tmp = g1 .* g1;')
_cg_rep('    curlFpG(GR1, GR2);\n    cout << "relaxfix',
        '    { real[int] g1(nt), g2(nt); curlFpG(g1, g2); real[int] q2 = g2 .* g2; tmp = g1 .* g1; tmp += q2; for [i, v : tmp] rhoG[i] = sqrt(v)/b; }   // rho_G of the relaxed state\n    cout << "relaxfix')
_cg_rep('    curlFpG(GR1, GR2);   // incompatibility present at the end of relaxation is the reference',
        '    { real[int] g1(nt), g2(nt); curlFpG(g1, g2); real[int] q2 = g2 .* g2; tmp = g1 .* g1; tmp += q2; for [i, v : tmp] rhoG[i] = sqrt(v)/b; }   // rho_G of the relaxed state (no reference)')
_cg_rep('// GND: G = Fp Curl Fp, rho_G = |G - G_R|/b, G_R = G at the end of the relaxation (<checkpoint>.gnd)',
        '// GND (Cermelli-Gurtin, planar): g = c/Jp, c_i = d1 Fp_i2 - d2 Fp_i1, rho_G = |g|/b at all times (no reference; the .gnd side file holds zeros)')
_cg_rep('// lumped P1 mass; element mean of c; G_i = Fp_ij c_j (reference curl, Cermelli-Gurtin G = Fp Curl Fp)',
        '// lumped P1 mass; element mean of c; Cermelli-Gurtin G = Jp^-1 Fp (Curl Fp)^T = e3 (x) c/Jp in plane strain: g = c/Jp')
_cg_rep('  // GND: |Fp Curl Fp - G_R|/b', '  // GND: rho_G = |c/Jp|/b (Cermelli-Gurtin)')
_cg_rep('thR = th; rhoG = 0;   // GND counts lattice curvature produced after relaxation only', 'thR = th;')
_cg_rep('\\"gnd\\": \\"curl Fp\\"', '\\"gnd\\": \\"Cermelli-Gurtin |c/Jp|/b, no reference\\"')
open(dst, "w").write(s)
print("wrote", dst)
