// dual.hpp -- minimal forward-mode automatic differentiation for the element kernels of kwc3d.cpp.
// D<N, S>: value v and N directional derivatives d[0..N-1] of scalar type S (S = double, or D<M> for second derivatives).
// Every operation applies its chain rule, so derivatives are those of the code as written, exact up to round-off.
// Checked against autodiff and CoDiPack (agreement 1e-14) and finite differences (checks/ad_check).
#pragma once
#include <cmath>

template <int N, class S = double> struct D {
  S v; S d[N];
  D() {}
  D(double x) : v(x) { for (int i = 0; i < N; i++) d[i] = 0.; }
  template <int M, class Q> D(const D<M, Q> &x) = delete;
};
template <class S> inline double valOf(const S &x) { return x; }
template <int N, class S> inline double valOf(const D<N, S> &x) { return valOf(x.v); }

#define K3_D template <int N, class S> inline D<N, S>
K3_D operator+(const D<N, S> &a, const D<N, S> &b) { D<N, S> r; r.v = a.v + b.v; for (int i = 0; i < N; i++) r.d[i] = a.d[i] + b.d[i]; return r; }
K3_D operator-(const D<N, S> &a, const D<N, S> &b) { D<N, S> r; r.v = a.v - b.v; for (int i = 0; i < N; i++) r.d[i] = a.d[i] - b.d[i]; return r; }
K3_D operator*(const D<N, S> &a, const D<N, S> &b) { D<N, S> r; r.v = a.v*b.v; for (int i = 0; i < N; i++) r.d[i] = a.d[i]*b.v + a.v*b.d[i]; return r; }
K3_D operator/(const D<N, S> &a, const D<N, S> &b) { D<N, S> r; S ib = 1./b.v; r.v = a.v*ib; for (int i = 0; i < N; i++) r.d[i] = (a.d[i] - r.v*b.d[i])*ib; return r; }
K3_D operator-(const D<N, S> &a) { D<N, S> r; r.v = -a.v; for (int i = 0; i < N; i++) r.d[i] = -a.d[i]; return r; }
K3_D operator+(const D<N, S> &a, double b) { D<N, S> r = a; r.v = r.v + b; return r; }
K3_D operator+(double b, const D<N, S> &a) { return a + b; }
K3_D operator-(const D<N, S> &a, double b) { D<N, S> r = a; r.v = r.v - b; return r; }
K3_D operator-(double b, const D<N, S> &a) { D<N, S> r = -a; r.v = r.v + b; return r; }
K3_D operator*(const D<N, S> &a, double b) { D<N, S> r; r.v = a.v*b; for (int i = 0; i < N; i++) r.d[i] = a.d[i]*b; return r; }
K3_D operator*(double b, const D<N, S> &a) { return a*b; }
K3_D operator/(const D<N, S> &a, double b) { return a*(1./b); }
K3_D operator/(double a, const D<N, S> &b) { D<N, S> r; S ib = 1./b.v; r.v = a*ib; for (int i = 0; i < N; i++) r.d[i] = -(r.v*b.d[i])*ib; return r; }
template <int N, class S> inline D<N, S> &operator+=(D<N, S> &a, const D<N, S> &b) { a = a + b; return a; }
template <int N, class S> inline D<N, S> &operator-=(D<N, S> &a, const D<N, S> &b) { a = a - b; return a; }
template <int N, class S> inline D<N, S> &operator+=(D<N, S> &a, double b) { a.v = a.v + b; return a; }
// f(a) with f = fv, f' = df (both of type S)
K3_D chainD(const D<N, S> &a, const S &fv, const S &df) { D<N, S> r; r.v = fv; for (int i = 0; i < N; i++) r.d[i] = df*a.d[i]; return r; }
using std::sin; using std::cos; using std::exp; using std::log; using std::sqrt;
K3_D sin(const D<N, S> &a) { return chainD(a, sin(a.v), cos(a.v)); }
K3_D cos(const D<N, S> &a) { return chainD(a, cos(a.v), -sin(a.v)); }
K3_D exp(const D<N, S> &a) { S e = exp(a.v); return chainD(a, e, e); }
K3_D log(const D<N, S> &a) { return chainD(a, log(a.v), 1./a.v); }
K3_D sqrt(const D<N, S> &a) { S s = sqrt(a.v); return chainD(a, s, 0.5/s); }
#undef K3_D
