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
#include "ff++.hpp"
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

template <class T> void cs(const T &q2, T &c, T &s) {
  double v = valOf(q2);
  if (fabs(v) < 1.) {   // series: c = sum q2^k/(2k)!, s = sum q2^k/(2k+1)!
    T tc = T(1.), ts = T(1.); c = T(1.); s = T(1.);
    for (int k = 1; k <= 12; k++) { tc = tc*q2/double((2*k - 1)*(2*k)); ts = ts*q2/double((2*k)*(2*k + 1)); c = c + tc; s = s + ts; }
  } else if (v > 0) { T r = sqrt(q2); T e = exp(r), ei = 1./e; c = 0.5*(e + ei); s = 0.5*(e - ei)/r; }
  else { T r = sqrt(-q2); c = cos(r); s = sin(r)/r; }
}
// out: PF11 PF12 PF21 PF22 psth ta tb tc tstar (and psi in out[9])
template <class T> void emap2(int k, const T *in, const double *fp, T *out) {
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
  T P11 = u11*S11 + u12*S12, P12 = u11*S12 + u12*S22 - ska, P21 = u21*S11 + u22*S12 + ska, P22 = u21*S12 + u22*S22;
  T b11 = ct*P11 - st*P21, b12 = ct*P12 - st*P22, b21 = st*P11 + ct*P21, b22 = st*P12 + ct*P22;
  out[0] = b11*fi11 + b12*fi12; out[1] = b11*fi21 + b12*fi22; out[2] = b21*fi11 + b22*fi12; out[3] = b21*fi21 + b22*fi22;
  out[4] = P11*u21 + P12*u22 - P21*u11 - P22*u12;
  T M11 = u11*P11 + u21*P21, M12 = u11*P12 + u21*P22, M21 = u12*P11 + u22*P21, M22 = u12*P12 + u22*P22;
  for (int a = 0; a < 3; a++) out[5 + a] = M11*Aa[4*a] + M12*Aa[4*a+1] + M21*Aa[4*a+2] + M22*Aa[4*a+3];
  out[8] = M21 - M12;
  out[9] = psi;
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
    typedef D<9> T; T in[9], o[10]; double fp[4];
    for (int j = 0; j < 9; j++) { in[j] = T(I(k, j)); in[j].d[j] = 1.; }
    for (int c = 0; c < 4; c++) fp[c] = I(k, 9 + c);
    emap2<T>(k, in, fp, o);
    for (int i = 0; i < 5; i++) for (int j = 0; j < 9; j++) S(k, i*9 + j) = o[i].d[j];
    for (int i = 5; i < 9; i++) for (int j = 0; j < 9; j++) Q(k, (i - 5)*9 + j) = o[i].d[j];
  }
  return 0;
}
// values only (checks): out nt x 10 (the 9 outputs and psi)
double k2dMap(KNM<double> *const &In, KNM<double> *const &Out) {
  KNM<double> &I = *In, &O = *Out;
  for (int k = 0; k < nt; k++) {
    double in[9], fp[4], o[10];
    for (int j = 0; j < 9; j++) in[j] = I(k, j); for (int c = 0; c < 4; c++) fp[c] = I(k, 9 + c);
    emap2<double>(k, in, fp, o);
    for (int c = 0; c < 10; c++) O(k, c) = o[c];
  }
  return 0;
}
double k2dThreads(const long &n) {
#ifdef _OPENMP
  if (n > 0) omp_set_num_threads(n); return omp_get_max_threads();
#else
  return 1;
#endif
}
static void Load_Init() {
  Global.Add("k2dSetup", "(", new OneOperator2_<double, KNM<double>*, KN<double>*>(k2dSetup));
  Global.Add("k2dJac", "(", new OneOperator3_<double, KNM<double>*, KNM<double>*, KNM<double>*>(k2dJac));
  Global.Add("k2dMap", "(", new OneOperator2_<double, KNM<double>*, KNM<double>*>(k2dMap));
  Global.Add("k2dThreads", "(", new OneOperator1_<double, long>(k2dThreads));
}
LOADFUNC(Load_Init)
