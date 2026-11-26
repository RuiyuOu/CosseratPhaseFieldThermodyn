// kwc2d.cpp -- exact derivatives of the 2D element map of cosserat_kwc.edp (FreeFem plugin), by forward-mode AD (dual.hpp).
// Build: ff-c++ kwc2d.cpp "CXXFLAGS=-fopenmp -O3" -lgomp
//
// Element map (same formulas as trialFp / finiteStress / eosStress / taus of the solver):
//   inputs  k = (F11-1, F12, F21, F22-1, theta) at a point, q = (dgamma_1..3, dw*), committed Fp (2x2)
//   outputs (dpsi/dF11, dpsi/dF12, dpsi/dF21, dpsi/dF22, dpsi/dtheta, tau_1, tau_2, tau_3, tau*)
// k2dJac returns d(outputs)/d(k, q) in the solver's layout: DSF(:, i*9 + j) for outputs i = 0..4, DTQ(:, (i-5)*9 + j) for i = 5..8,
// inputs j = 0..8 = (e11, hy1, hx2, e22, thc, dpa, dpb, dpc, dwst).  exp of the traceless 2x2 Lp: exp(L) = c(q2) I + s(q2) L with
// c = cosh(sqrt q2), s = sinh(sqrt q2)/sqrt q2 written as power series in q2 for |q2| < 1 (entire in q2: smooth through q2 = 0 for AD)
// and in closed form (cosh/sinh or cos/sin) otherwise.
#include <algorithm>
#include <functional>
#include <vector>
#include "ff++.hpp"
#include "AFunction_ext.hpp"
#include "dual.hpp"
#include <vector>
#include <cmath>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace std;

namespace k2 {
int nt = 0;
vector<double> CS;   // nt x 21: C11 C22 C12 C16 C26 C66 B11 B22 B12, A11a A12a A21a A22a, A11b .. A22b, A11c .. A22c
double muc = 75.4, K0 = 137.07, Kp = 5., Clxx = 220.3;
vector<double> LAM;  // v2 solver only (k2dLam): Cosserat multiplier per element, psi += lambda a; empty = none

template <class T> void cs(const T &q2, T &c, T &s) {
  double v = valOf(q2);
  if (fabs(v) < 1.) {   // series: c = sum q2^k/(2k)!, s = sum q2^k/(2k+1)!
    T tc = T(1.), ts = T(1.); c = T(1.); s = T(1.);
    for (int k = 1; k <= 12; k++) { tc = tc*q2/double((2*k - 1)*(2*k)); ts = ts*q2/double((2*k)*(2*k + 1)); c = c + tc; s = s + ts; }
  } else if (v > 0) { T r = sqrt(q2); T e = exp(r), ei = 1./e; c = 0.5*(e + ei); s = 0.5*(e - ei)/r; }
  else { T r = sqrt(-q2); c = cos(r); s = sin(r)/r; }
}
// out: PF11 PF12 PF21 PF22 psth ta tb tc tstar (psi in out[9], the Cosserat mismatch a in out[10]).
// lam: differentiable multiplier (k2dJacL), or null -> the stored LAM[k] as a constant
template <class T> void emap2(int k, const T *in, const double *fp, T *out, const T *lam = 0) {
  const double *e = &CS[21*k];
  const double C11 = e[0], C22 = e[1], C12 = e[2], C16 = e[3], C26 = e[4], C66 = e[5], B11 = e[6], B22 = e[7], B12 = e[8];
  const double *Aa = e + 9;   // A11,A12,A21,A22 for system a at Aa + 4a
  T dp[3] = {in[5], in[6], in[7]}, dw = in[8];
  T L11 = 0., L12 = 0., L21 = 0., L22 = 0.;
  for (int a = 0; a < 3; a++) { L11 += Aa[4*a]*dp[a]; L12 += Aa[4*a+1]*dp[a]; L21 += Aa[4*a+2]*dp[a]; L22 += Aa[4*a+3]*dp[a]; }
  L12 = L12 - dw; L21 += dw;
  T q2 = L11*L11 + L12*L21, c, s; cs(q2, c, s);
  T X11 = c + s*L11, X12 = s*L12, X21 = s*L21, X22 = c + s*L22;
  T tf11 = X11*fp[0] + X12*fp[2], tf12 = X11*fp[1] + X12*fp[3], tf21 = X21*fp[0] + X22*fp[2], tf22 = X21*fp[1] + X22*fp[3];
  T dtf = tf11*tf22 - tf12*tf21;
  T fi11 = tf22/dtf, fi22 = tf11/dtf, fi12 = -tf12/dtf, fi21 = -tf21/dtf;
  T F11 = in[0] + 1., F12 = in[1], F21 = in[2], F22 = in[3] + 1.;
  T A11 = F11*fi11 + F12*fi21, A12 = F11*fi12 + F12*fi22, A21 = F21*fi11 + F22*fi21, A22 = F21*fi12 + F22*fi22;
  T ct = cos(in[4]), st = sin(in[4]);
  T u11 = ct*A11 + st*A21, u12 = ct*A12 + st*A22, u21 = ct*A21 - st*A11, u22 = ct*A22 - st*A12;
  T Eb11 = 0.5*(u11*u11 + u21*u21 - 1.), Eb22 = 0.5*(u12*u12 + u22*u22 - 1.), Eb12 = 0.5*(u11*u12 + u21*u22);
  T J = u11*u22 - u12*u21;
  T Jc = (valOf(J) > 1e-6) ? J : T(1e-6);
  T lJ = log(Jc), Dd = exp(-2./3*lJ);
  T Cb11 = 2.*Eb11 + 1., Cb22 = 2.*Eb22 + 1., Cb12 = 2.*Eb12;
  T Et11 = 0.5*(Dd*Cb11 - 1.), Et22 = 0.5*(Dd*Cb22 - 1.), Et12 = 0.5*Dd*Cb12, Et33 = 0.5*(Dd - 1.);
  T q11 = C11*Et11 + C12*Et22 + 2*C16*Et12, q22 = C12*Et11 + C22*Et22 + 2*C26*Et12, q12 = C16*Et11 + C26*Et22 + 2*C66*Et12;
  T bE = B11*Et11 + B22*Et22 + 2*B12*Et12;
  T psi = 0.5*(q11*Et11 + q22*Et22 + 2.*q12*Et12) + Et33*bE + 0.5*Clxx*Et33*Et33;
  T T11 = q11 + B11*Et33, T22 = q22 + B22*Et33, T12 = q12 + B12*Et33, T33 = Clxx*Et33 + bE;
  double eta = 1.5*(Kp - 1);
  T x = exp(lJ/3.), omx = 1. - x, ee = exp(eta*omx);
  T Peos = 3*K0*omx/(x*x)*ee;
  psi += 9*K0/(eta*eta)*(1. - (1. - eta*omx)*ee);
  T J2 = Jc*Jc, Ci11 = Cb22/J2, Ci22 = Cb11/J2, Ci12 = -(Cb12/J2);
  T trT = (T11*Cb11 + T22*Cb22 + 2.*T12*Cb12 + T33)/3.;
  T PJ = Peos*Jc;
  T S11 = Dd*(T11 - trT*Ci11) - PJ*Ci11, S22 = Dd*(T22 - trT*Ci22) - PJ*Ci22, S12 = Dd*(T12 - trT*Ci12) - PJ*Ci12;
  T aa = 0.5*(u21 - u12); psi += 2*muc*aa*aa; T ska = 2*muc*aa;
  if (lam) { psi += (*lam)*aa; ska = ska + 0.5*(*lam); }
  else if (!LAM.empty() && LAM[k] != 0.) { psi += LAM[k]*aa; ska = ska + 0.5*LAM[k]; }
  T P11 = u11*S11 + u12*S12, P12 = u11*S12 + u12*S22 - ska, P21 = u21*S11 + u22*S12 + ska, P22 = u21*S12 + u22*S22;
  T b11 = ct*P11 - st*P21, b12 = ct*P12 - st*P22, b21 = st*P11 + ct*P21, b22 = st*P12 + ct*P22;
  out[0] = b11*fi11 + b12*fi12; out[1] = b11*fi21 + b12*fi22; out[2] = b21*fi11 + b22*fi12; out[3] = b21*fi21 + b22*fi22;
  out[4] = P11*u21 + P12*u22 - P21*u11 - P22*u12;
  T M11 = u11*P11 + u21*P21, M12 = u11*P12 + u21*P22, M21 = u12*P11 + u22*P21, M22 = u12*P12 + u22*P22;
  for (int a = 0; a < 3; a++) out[5 + a] = M11*Aa[4*a] + M12*Aa[4*a+1] + M21*Aa[4*a+2] + M22*Aa[4*a+3];
  out[8] = M21 - M12;
  out[9] = psi; out[10] = aa;
}
} // namespace k2
using namespace k2;

// set the per-element constants (nt x 21) and the parameters (muc, K0, K', C_zzzz)
double k2dSetup(KNM<double> *const &C, KN<double> *const &p) {
  KNM<double> &A = *C; nt = A.N(); CS.resize(21*nt);
  for (int k = 0; k < nt; k++) for (int c = 0; c < 21; c++) CS[21*k + c] = A(k, c);
  muc = (*p)[0]; K0 = (*p)[1]; Kp = (*p)[2]; Clxx = (*p)[3];
  return nt;
}
// in: nt x 13 (e11 hy1 hx2 e22 thc dpa dpb dpc dwst fp11 fp12 fp21 fp22) -> DSF (nt x 45), DTQ (nt x 36)
double k2dJac(KNM<double> *const &In, KNM<double> *const &DSF, KNM<double> *const &DTQ) {
  KNM<double> &I = *In, &S = *DSF, &Q = *DTQ;
#pragma omp parallel for schedule(static)
  for (int k = 0; k < nt; k++) {
    typedef D<9> T; T in[9], o[11]; double fp[4];
    for (int j = 0; j < 9; j++) { in[j] = T(I(k, j)); in[j].d[j] = 1.; }
    for (int c = 0; c < 4; c++) fp[c] = I(k, 9 + c);
    emap2<T>(k, in, fp, o);
    for (int i = 0; i < 5; i++) for (int j = 0; j < 9; j++) S(k, i*9 + j) = o[i].d[j];
    for (int i = 5; i < 9; i++) for (int j = 0; j < 9; j++) Q(k, (i - 5)*9 + j) = o[i].d[j];
  }
  return 0;
}
// v2 solver, monolithic multiplier (-almc 4): as k2dJac with the element multiplier as a 10th differentiated input.
// DAL (nt x 18): [0..8] da/d(k, q), [9..13] d(stress outputs 0..4)/dlambda, [14..17] d(tau_1..3, tau*)/dlambda
double k2dJacL(KNM<double> *const &In, KNM<double> *const &DSF, KNM<double> *const &DTQ, KNM<double> *const &DAL) {
  KNM<double> &I = *In, &S = *DSF, &Q = *DTQ, &A = *DAL;
#pragma omp parallel for schedule(static)
  for (int k = 0; k < nt; k++) {
    typedef D<10> T; T in[9], o[11], lam; double fp[4];
    for (int j = 0; j < 9; j++) { in[j] = T(I(k, j)); in[j].d[j] = 1.; }
    lam = T(LAM.empty() ? 0. : LAM[k]); lam.d[9] = 1.;
    for (int c = 0; c < 4; c++) fp[c] = I(k, 9 + c);
    emap2<T>(k, in, fp, o, &lam);
    for (int i = 0; i < 5; i++) for (int j = 0; j < 9; j++) S(k, i*9 + j) = o[i].d[j];
    for (int i = 5; i < 9; i++) for (int j = 0; j < 9; j++) Q(k, (i - 5)*9 + j) = o[i].d[j];
    for (int j = 0; j < 9; j++) A(k, j) = o[10].d[j];
    for (int i = 0; i < 9; i++) A(k, 9 + i) = o[i].d[9];
  }
  return 0;
}
// ==== direct assembly of the tangent matrix (v2 solver, -asm 1) =====================================================
// Fixed mesh and reference configuration: the sparsity pattern, the kinematic operators G_q (Gauss points) and G_c (centroid), the
// element -> matrix-entry map and the constant part (mass, viscosity, boundary terms) are built once.  Per Newton iteration only the
// Gauss-point blocks change:  K_e = (A_e/3) sum_q G_q^T (HQ_q G_q + CQ_q G_c) + A_e gradN^T H_kwc gradN.
namespace k2a {
  int nel = 0; long nnz = 0, ndof = 0;
  std::vector<int> conn;            // nel x 15 global dofs (u1: 6, u2: 6, theta: 3)
  std::vector<double> geo;          // nel x 55: gbX (24), gbY (24), P1 gradients gx (3), gy (3), area
  std::vector<int> pos;             // nel x 225: index of entry (a, b) in the coefficient array of the FreeFEM matrix
  std::vector<long> mi, mj;         // the matrix entries in that order (probe)
  std::vector<double> base, jg;     // constant part; Gauss-point Hessian part of the current iterate
  std::vector<double> pendH, pendC; long pendN = 0;   // Gauss-point blocks handed over before the set-up
  std::vector<long> keyS; std::vector<int> keyI;   // sorted (row*ndof + column) keys of the pattern and their entry index (ASMCHECK)
  const double gT[4][3] = {{.5, .5, 0.}, {0., .5, .5}, {.5, 0., .5}, {1./3, 1./3, 1./3}};
  inline void opG(const double *g, int q, double G[5][15]) {   // rows: du1/dx, du1/dy, du2/dx, du2/dy, theta
    for (int i = 0; i < 5; i++) for (int a = 0; a < 15; a++) G[i][a] = 0.;
    for (int a = 0; a < 6; a++) { G[0][a] = g[q*6 + a]; G[1][a] = g[24 + q*6 + a]; G[2][6 + a] = g[q*6 + a]; G[3][6 + a] = g[24 + q*6 + a]; }
    for (int a = 0; a < 3; a++) G[4][12 + a] = gT[q][a];
  }
}
static void asmJcore(const std::function<double(int,int)> &H, const std::function<double(int,int)> &Cq);
// unique (row, column) pairs of the element couplings
double k2dAsmPattern(KNM<long> *const &C, KN<long> *const &I, KN<long> *const &J) {
  int ne = C->N(); std::vector<std::pair<long,long> > v; v.reserve((size_t)ne*225);
  for (int k = 0; k < ne; k++) for (int a = 0; a < 15; a++) for (int b = 0; b < 15; b++) v.push_back(std::make_pair((*C)(k, a), (*C)(k, b)));
  std::sort(v.begin(), v.end()); v.erase(std::unique(v.begin(), v.end()), v.end());
  I->resize(v.size()); J->resize(v.size());
  for (long e = 0; e < (long)v.size(); e++) { (*I)[e] = v[e].first; (*J)[e] = v[e].second; }
  return (double)v.size();
}
// conn, geo, entry order (I, J) of the FreeFEM matrix, constant matrix as triplets (Ik, Jk, Ck).  Returns the number of constant
// entries outside the pattern (must be 0).
class K2AsmSetup_Op : public E_F0mps { public:
  Expression e[7];
  K2AsmSetup_Op(const basicAC_F0 &args) { args.SetNameParam(0, 0, 0); for (int i = 0; i < 7; i++) e[i] = (i < 2) ? (i == 0 ? to< KNM<long>* >(args[i]) : to< KNM<double>* >(args[i])) : ((i == 6) ? to< KN<double>* >(args[i]) : to< KN<long>* >(args[i])); }
  static ArrayOfaType typeargs() { return ArrayOfaType(atype< KNM<long>* >(), atype< KNM<double>* >(), atype< KN<long>* >(), atype< KN<long>* >(), atype< KN<long>* >(), atype< KN<long>* >(), atype< KN<double>* >(), false); }
  static E_F0 *f(const basicAC_F0 &args) { return new K2AsmSetup_Op(args); }
  AnyType operator()(Stack s) const {
    using namespace k2a;
    KNM<long> &C = *GetAny< KNM<long>* >((*e[0])(s)); KNM<double> &Gm = *GetAny< KNM<double>* >((*e[1])(s));
    KN<long> &I = *GetAny< KN<long>* >((*e[2])(s)), &J = *GetAny< KN<long>* >((*e[3])(s)), &Ik = *GetAny< KN<long>* >((*e[4])(s)), &Jk = *GetAny< KN<long>* >((*e[5])(s));
    KN<double> &Ck = *GetAny< KN<double>* >((*e[6])(s));
    nel = C.N(); nnz = I.N(); conn.resize((size_t)nel*15); geo.resize((size_t)nel*55); pos.resize((size_t)nel*225);
    for (int k = 0; k < nel; k++) { for (int a = 0; a < 15; a++) conn[k*15 + a] = (int)C(k, a); for (int c = 0; c < 55; c++) geo[(size_t)k*55 + c] = Gm(k, c); }
    ndof = 0; for (long e2 = 0; e2 < nnz; e2++) ndof = std::max(ndof, std::max(I[e2], J[e2]) + 1);
    std::vector<std::pair<std::pair<long,long>, int> > v(nnz);
    mi.resize(nnz); mj.resize(nnz);
    for (long e2 = 0; e2 < nnz; e2++) { v[e2] = std::make_pair(std::make_pair(I[e2], J[e2]), (int)e2); mi[e2] = I[e2]; mj[e2] = J[e2]; }
    std::sort(v.begin(), v.end());
    keyS.resize(nnz); keyI.resize(nnz);
    for (long e2 = 0; e2 < nnz; e2++) { keyS[e2] = v[e2].first.first*ndof + v[e2].first.second; keyI[e2] = v[e2].second; }
    long miss = 0;
    auto find = [&](long i, long j) -> int {
      auto it = std::lower_bound(v.begin(), v.end(), std::make_pair(std::make_pair(i, j), -1));
      if (it == v.end() || it->first.first != i || it->first.second != j) return -1;
      return it->second; };
    for (int k = 0; k < nel; k++) for (int a = 0; a < 15; a++) for (int b = 0; b < 15; b++) {
      int f = find(conn[k*15 + a], conn[k*15 + b]); if (f < 0) { miss++; f = 0; } pos[(size_t)k*225 + a*15 + b] = f; }
    base.assign(nnz, 0.); jg.assign(nnz, 0.);
    for (long e2 = 0; e2 < Ik.N(); e2++) { int f = find(Ik[e2], Jk[e2]); if (f < 0) { if (Ck[e2] != 0.) miss++; } else base[f] += Ck[e2]; }
    if (pendN == nel) {   // Gauss-point blocks handed over before the set-up
      asmJcore([&](int k, int c) { return pendH[(size_t)k*75 + c]; }, [&](int k, int c) { return pendC[(size_t)k*75 + c]; });
      pendH.clear(); pendC.clear(); pendN = 0;
    } else if (pendN != 0) miss += 1000000000;   // blocks of another size: refuse
    return SetAny<double>((double)miss);
  }
};
class K2AsmSetup : public OneOperator { public:
  K2AsmSetup() : OneOperator(atype<double>(), atype< KNM<long>* >(), atype< KNM<double>* >(), atype< KN<long>* >(), atype< KN<long>* >(), atype< KN<long>* >(), atype< KN<long>* >(), atype< KN<double>* >()) {}
  E_F0 *code(const basicAC_F0 &args) const { return new K2AsmSetup_Op(args); }
};
// Gauss-point part: HQ, CQ (nt x 75; column q*25 + i*5 + j).  Blocks that arrive before k2dAsmSetup (the first tangent of a run is
// built before the first matrix assembly, where the set-up is made) are kept and processed at the end of the set-up.
static void asmJcore(const std::function<double(int,int)> &H, const std::function<double(int,int)> &Cq) {
  using namespace k2a;
  std::fill(jg.begin(), jg.end(), 0.);
  for (int k = 0; k < nel; k++) {
    const double *g = &geo[(size_t)k*55]; double Gc[5][15], G[5][15], T[5][15], Ke[225]; const double w = g[54]/3.;
    opG(g, 3, Gc); for (int e = 0; e < 225; e++) Ke[e] = 0.;
    for (int q = 0; q < 3; q++) {
      opG(g, q, G);
      for (int i = 0; i < 5; i++) for (int b = 0; b < 15; b++) {
        double t = 0.; for (int j = 0; j < 5; j++) t += H(k, q*25 + i*5 + j)*G[j][b] + Cq(k, q*25 + i*5 + j)*Gc[j][b]; T[i][b] = t*w; }
      for (int a = 0; a < 15; a++) for (int i = 0; i < 5; i++) { const double ga = G[i][a]; if (ga != 0.) for (int b = 0; b < 15; b++) Ke[a*15 + b] += ga*T[i][b]; }
    }
    const int *p = &pos[(size_t)k*225]; for (int e = 0; e < 225; e++) jg[p[e]] += Ke[e];
  }
}
double k2dAsmJ(KNM<double> *const &HQ, KNM<double> *const &CQ) {
  using namespace k2a;
  if (nel == 0) {   // no set-up yet: keep a copy
    long n = HQ->N(); pendN = n; pendH.resize((size_t)n*75); pendC.resize((size_t)n*75);
    for (int k = 0; k < (int)n; k++) for (int c = 0; c < 75; c++) { pendH[(size_t)k*75 + c] = (*HQ)(k, c); pendC[(size_t)k*75 + c] = (*CQ)(k, c); }
    return 0;
  }
  asmJcore([&](int k, int c) { return (*HQ)(k, c); }, [&](int k, int c) { return (*CQ)(k, c); });
  return 0;
}
// full matrix values: constant part + (useJ ? Gauss-point part : 0) + KWC curvature tangent (HK nt x 3: H11, H12, H22; nt x 4: + Hd)
double k2dAsmA(KNM<double> *const &HK, long const &useJ, KN<double> *const &vals) {
  using namespace k2a;
  KN<double> &V = *vals; ffassert(V.N() == nnz);
  const bool hk4 = HK->M() >= 4;   // fourth column: antisymmetric part (H12 - H21)/2 of the KWC tangent (-kwctan 3)
  if (useJ) for (long e = 0; e < nnz; e++) V[e] = base[e] + jg[e]; else for (long e = 0; e < nnz; e++) V[e] = base[e];
  for (int k = 0; k < nel; k++) {
    const double *g = &geo[(size_t)k*55]; const double *gx = g + 48, *gy = g + 51, A = g[54];
    const double h11 = (*HK)(k, 0), h12 = (*HK)(k, 1), h22 = (*HK)(k, 2), hd = (hk4 ? (*HK)(k, 3) : 0.); const int *p = &pos[(size_t)k*225];
    for (int a = 0; a < 3; a++) for (int b = 0; b < 3; b++)   // row a: test function zeta_a, column b: theta_b
      V[p[(12 + a)*15 + 12 + b]] += A*(h11*gx[b]*gx[a] + h12*(gx[b]*gy[a] + gy[b]*gx[a]) + h22*gy[b]*gy[a] + hd*(gy[b]*gx[a] - gx[b]*gy[a]));
  }
  return 0;
}
// z^T A x with the stored entry order: compared with the FreeFEM product after every assembly (guards the entry order)
double k2dAsmProbe(KN<double> *const &vals, KN<double> *const &x, KN<double> *const &z) {
  using namespace k2a; double s = 0.;
  for (long e = 0; e < nnz; e++) s += (*z)[mi[e]]*(*vals)[e]*(*x)[mj[e]];
  return s;
}
// sum of |z_i A_ij x_j|: the scale for the comparison (z'Ax itself can pass through zero)
double k2dAsmProbeAbs(KN<double> *const &vals, KN<double> *const &x, KN<double> *const &z) {
  using namespace k2a; double s = 0.;
  for (long e = 0; e < nnz; e++) s += fabs((*z)[mi[e]]*(*vals)[e]*(*x)[mj[e]]);
  return s;
}
// ASMCHECK: compare a matrix given as triplets (I, J, C) (the form assembly) with the direct values V on the pattern.
// Returns max |difference| / max |V| over the rows without a Dirichlet diagonal (|diag| < 1e20); an entry of (I, J, C) outside the
// pattern counts with its full value, an entry of the pattern missing from (I, J, C) with the full value of V.
double k2dAsmCmp(KN<long> *const &I, KN<long> *const &J, KN<double> *const &C, KN<double> *const &vals) {
  using namespace k2a; KN<double> &V = *vals;
  std::vector<char> dir(ndof, 0), seen(nnz, 0);
  for (long e = 0; e < nnz; e++) if (mi[e] == mj[e] && fabs(V[e]) >= 1e20) dir[mi[e]] = 1;
  double vmax = 0., dmax = 0.;
  for (long e = 0; e < nnz; e++) if (!dir[mi[e]]) vmax = std::max(vmax, fabs(V[e]));
  for (long t = 0; t < I->N(); t++) {
    long i = (*I)[t], j = (*J)[t]; if (i < 0 || i >= ndof || j < 0 || j >= ndof) { dmax = 1e300; continue; }
    if (dir[i]) continue;
    long key = i*ndof + j; auto it = std::lower_bound(keyS.begin(), keyS.end(), key);
    if (it == keyS.end() || *it != key) { dmax = std::max(dmax, fabs((*C)[t])); continue; }
    int e = keyI[it - keyS.begin()]; seen[e] = 1; dmax = std::max(dmax, fabs((*C)[t] - V[e]));
  }
  for (long e = 0; e < nnz; e++) if (!seen[e] && !dir[mi[e]]) dmax = std::max(dmax, fabs(V[e]));
  return dmax/std::max(vmax, 1e-300);
}
// values only (checks): out nt x 10 (the 9 outputs and psi), or nt x 11 (+ a)
double k2dMap(KNM<double> *const &In, KNM<double> *const &Out) {
  KNM<double> &I = *In, &O = *Out;
  for (int k = 0; k < nt; k++) {
    double in[9], fp[4], o[11];
    for (int j = 0; j < 9; j++) in[j] = I(k, j); for (int c = 0; c < 4; c++) fp[c] = I(k, 9 + c);
    emap2<double>(k, in, fp, o);
    for (int c = 0; c < O.M() && c < 11; c++) O(k, c) = o[c];
  }
  return 0;
}
// v2 solver (-almc): per-element Cosserat multiplier lambda (the base solver never calls this)
double k2dLam(KN<double> *const &L) { LAM.resize(L->N()); for (long k = 0; k < L->N(); k++) LAM[k] = (*L)[k]; return 0; }
double k2dThreads(const long &n) {
#ifdef _OPENMP
  if (n > 0) omp_set_num_threads(n); return omp_get_max_threads();
#else
  return 1;
#endif
}
// ---------------- phase field (fast assembly, r8): the forms vphA, vphR, vph2, vgph, vDph of the solver on a fixed pattern ----------------
// Same 7-point rule as FreeFEM's int2d default (qf5pT). P1 phi, phik, phin, T; P0 N, rho_S.
namespace kph {
  int nel = 0; long nnz = 0, ndof = 0;
  std::vector<int> conn, pos; std::vector<double> area, gx, gy;
  double L[7][3], W[7];
  void rule() {
    const double s = sqrt(15.), a1 = (6 - s)/21, b1 = (9 + 2*s)/21, a2 = (6 + s)/21, b2 = (9 - 2*s)/21;
    const double P[7][3] = {{1./3, 1./3, 1./3}, {a1, a1, b1}, {a1, b1, a1}, {b1, a1, a1}, {a2, a2, b2}, {a2, b2, a2}, {b2, a2, a2}};
    for (int q = 0; q < 7; q++) { for (int i = 0; i < 3; i++) L[q][i] = P[q][i]; W[q] = (q == 0) ? 9./40 : (q < 4 ? (155 - s)/1200 : (155 + s)/1200); }
  }
}
// set-up: C (nt x 3 P1 dofs), XY (nt x 6 vertex coordinates), I, J = entry order of the phase matrix ([I, J, C] = A). Returns the entries missed (0).
double k2dPhSetup(KNM<long> *const &C, KNM<double> *const &XY, KN<long> *const &I, KN<long> *const &J) {
  using namespace kph; rule();
  nel = C->N(); nnz = I->N(); conn.resize((size_t)nel*3); area.resize(nel); gx.resize((size_t)nel*3); gy.resize((size_t)nel*3); pos.resize((size_t)nel*9);
  ndof = 0; for (long e = 0; e < nnz; e++) ndof = std::max(ndof, std::max((*I)[e], (*J)[e]) + 1);
  for (int k = 0; k < nel; k++) ndof = std::max(ndof, std::max((*C)(k, 0), std::max((*C)(k, 1), (*C)(k, 2))) + 1);
  std::vector<std::pair<long,long> > key(nnz); for (long e = 0; e < nnz; e++) key[e] = std::make_pair((*I)[e]*ndof + (*J)[e], e);
  std::sort(key.begin(), key.end());
  long miss = 0;
  for (int k = 0; k < nel; k++) {
    for (int a = 0; a < 3; a++) conn[k*3 + a] = (int)(*C)(k, a);
    const double x0 = (*XY)(k,0), y0 = (*XY)(k,1), x1 = (*XY)(k,2), y1 = (*XY)(k,3), x2 = (*XY)(k,4), y2 = (*XY)(k,5);
    const double d = (x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0); area[k] = 0.5*fabs(d);
    gx[k*3+0] = (y1 - y2)/d; gx[k*3+1] = (y2 - y0)/d; gx[k*3+2] = (y0 - y1)/d;
    gy[k*3+0] = (x2 - x1)/d; gy[k*3+1] = (x0 - x2)/d; gy[k*3+2] = (x1 - x0)/d;
    for (int a = 0; a < 3; a++) for (int b = 0; b < 3; b++) {
      const long kk = (long)conn[k*3+a]*ndof + conn[k*3+b];
      auto it = std::lower_bound(key.begin(), key.end(), std::make_pair(kk, -1L));
      if (it == key.end() || it->first != kk) { miss++; pos[(size_t)k*9 + a*3 + b] = 0; } else pos[(size_t)k*9 + a*3 + b] = (int)it->second;
    }
  }
  return (double)miss;
}
// matrix values V (vphA, linearised at phik) and residual R (vphR, at phi).  F (ndof x 4): phi, phik, phin, T.  E (nt x 2): N, rho_S.
// P: etaPhi, dt, f0, sK, aK, T0, Q/R, Arrhenius floor, c_s = chi mu b^2/2.
double k2dPhAsm(KNM<double> *const &F, KNM<double> *const &E, KN<double> *const &P, KN<double> *const &V, KN<double> *const &R) {
  using namespace kph;
  const KNM<double> &f = *F; const KNM<double> &el = *E; const KN<double> &p = *P;
  const double eta = p[0], dt = p[1], f0 = p[2], sK = p[3], aK = p[4], T0 = p[5], QR = p[6], amin = p[7], cs = p[8], kK = f0*aK*aK;
  KN<double> &v = *V, &r = *R; ffassert(nel > 0 && v.N() == nnz && r.N() == ndof && el.N() == nel);
  for (long e = 0; e < nnz; e++) v[e] = 0.; for (long i = 0; i < ndof; i++) r[i] = 0.;
  for (int k = 0; k < nel; k++) {
    const int *c = &conn[(size_t)k*3]; const double A = area[k], N = el(k, 0), rho = el(k, 1);
    double ph[3], pk[3], pn[3], T[3]; for (int a = 0; a < 3; a++) { ph[a] = f(c[a], 0); pk[a] = f(c[a], 1); pn[a] = f(c[a], 2); T[a] = f(c[a], 3); }
    double Me[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0}, re[3] = {0, 0, 0};
    for (int q = 0; q < 7; q++) {
      const double *l = L[q]; double phq = 0, pkq = 0, pnq = 0, Tq = 0;
      for (int a = 0; a < 3; a++) { phq += l[a]*ph[a]; pkq += l[a]*pk[a]; pnq += l[a]*pn[a]; Tq += l[a]*T[a]; }
      const double arr = std::max(exp(QR*(1./Tq - 1./T0)), amin)*eta/dt;
      const double om = 1. - pkq, cq = arr + f0 + f0*sK*N*2./(om*om);
      const double fq = arr*(phq - pnq) + f0*phq + f0*sK*N*2.*phq/(1. - phq) - (f0 - cs*rho);
      const double w = W[q]*A;
      for (int a = 0; a < 3; a++) { re[a] += w*fq*l[a]; for (int b = 0; b < 3; b++) Me[a*3 + b] += w*cq*l[a]*l[b]; }
    }
    const double *Gx = &gx[(size_t)k*3], *Gy = &gy[(size_t)k*3];
    for (int a = 0; a < 3; a++) {
      double kr = 0;
      for (int b = 0; b < 3; b++) { const double ke = kK*A*(Gx[a]*Gx[b] + Gy[a]*Gy[b]); v[pos[(size_t)k*9 + a*3 + b]] += Me[a*3 + b] + ke; kr += ke*ph[b]; }
      r[c[a]] += re[a] + kr;
    }
  }
  return 0;
}
// element integrals divided by the area (vph2, vgph, vDph of the solver): O (nt x 3) = mean phi^2, mean g(phi), mean eta_phi(T)(phi - phin)^2/dt.
// F (ndof x 3): phi, phin, T.  P as for k2dPhAsm (etaPhi, dt, ..., Q/R, floor).
double k2dPhElem(KNM<double> *const &F, KN<double> *const &P, KNM<double> *const &O) {
  using namespace kph;
  const KNM<double> &f = *F; const KN<double> &p = *P; KNM<double> &o = *O; ffassert(nel > 0 && o.N() == nel && o.M() >= 3);
  const double eta = p[0], dt = p[1], T0 = p[5], QR = p[6], amin = p[7];
  for (int k = 0; k < nel; k++) {
    const int *c = &conn[(size_t)k*3]; double s2 = 0, sg = 0, sd = 0;
    for (int q = 0; q < 7; q++) {
      double phq = 0, pnq = 0, Tq = 0; for (int a = 0; a < 3; a++) { phq += L[q][a]*f(c[a], 0); pnq += L[q][a]*f(c[a], 1); Tq += L[q][a]*f(c[a], 2); }
      const double arr = std::max(exp(QR*(1./Tq - 1./T0)), amin), dp = phq - pnq;
      s2 += W[q]*phq*phq; sg += W[q]*(-2.*(phq + log(1. - phq))); sd += W[q]*eta*arr*dp*dp/dt;
    }
    o(k, 0) = s2; o(k, 1) = sg; o(k, 2) = sd;   // the weights sum to 1: these are the element means
  }
  return 0;
}
static void Load_Init() {
  Global.Add("k2dSetup", "(", new OneOperator2_<double, KNM<double>*, KN<double>*>(k2dSetup));
  Global.Add("k2dJac", "(", new OneOperator3_<double, KNM<double>*, KNM<double>*, KNM<double>*>(k2dJac));
  Global.Add("k2dJacL", "(", new OneOperator4_<double, KNM<double>*, KNM<double>*, KNM<double>*, KNM<double>*>(k2dJacL));
  Global.Add("k2dAsmPattern", "(", new OneOperator3_<double, KNM<long>*, KN<long>*, KN<long>*>(k2dAsmPattern));
  Global.Add("k2dAsmSetup", "(", new K2AsmSetup);
  Global.Add("k2dAsmJ", "(", new OneOperator2_<double, KNM<double>*, KNM<double>*>(k2dAsmJ));
  Global.Add("k2dAsmA", "(", new OneOperator3_<double, KNM<double>*, long, KN<double>*>(k2dAsmA));
  Global.Add("k2dAsmProbe", "(", new OneOperator3_<double, KN<double>*, KN<double>*, KN<double>*>(k2dAsmProbe));
  Global.Add("k2dAsmCmp", "(", new OneOperator4_<double, KN<long>*, KN<long>*, KN<double>*, KN<double>*>(k2dAsmCmp));
  Global.Add("k2dAsmProbeAbs", "(", new OneOperator3_<double, KN<double>*, KN<double>*, KN<double>*>(k2dAsmProbeAbs));
  Global.Add("k2dMap", "(", new OneOperator2_<double, KNM<double>*, KNM<double>*>(k2dMap));
  Global.Add("k2dLam", "(", new OneOperator1_<double, KN<double>*>(k2dLam));
  Global.Add("k2dPhSetup", "(", new OneOperator4_<double, KNM<long>*, KNM<double>*, KN<long>*, KN<long>*>(k2dPhSetup));
  Global.Add("k2dPhAsm", "(", new OneOperator5_<double, KNM<double>*, KNM<double>*, KN<double>*, KN<double>*, KN<double>*>(k2dPhAsm));
  Global.Add("k2dPhElem", "(", new OneOperator3_<double, KNM<double>*, KN<double>*, KNM<double>*>(k2dPhElem));
  Global.Add("k2dThreads", "(", new OneOperator1_<double, long>(k2dThreads));
}
LOADFUNC(Load_Init)
