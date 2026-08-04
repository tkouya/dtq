/*
 * src/qs_real.cpp
 *
 * Implementation of the qs_real (quad-single, 4 x float) class.
 * Constants and transcendentals delegate to the higher-precision qd_real
 * routines and the result is rounded back to ~96-bit (qs_real) precision.
 *
 * Created for the dtq package extension (2026).
 */
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>

#include "config.h"
#include <qd/qs_real.h>
#include <qd/qd_real.h>

using std::cout;
using std::cerr;
using std::endl;
using std::ostream;
using std::istream;
using std::string;

/*======== Helpers: qs_real <-> qd_real ========*/
/* Extract the top ~96 bits of `orig` into four IEEE-754 floats.
   The residual after each extraction is computed using qd_real arithmetic
   so that no precision is lost between successive 24-bit chunks
   (using a plain double residual would only keep 53 bits and would zero
   out the fourth float for values near 1.0).                           */
static qs_real qs_from_qd(const qd_real &orig) {
  qd_real rem = orig;
  float a0 = static_cast<float>(rem[0]);
  rem -= static_cast<double>(a0);
  float a1 = static_cast<float>(rem[0]);
  rem -= static_cast<double>(a1);
  float a2 = static_cast<float>(rem[0]);
  rem -= static_cast<double>(a2);
  float a3 = static_cast<float>(rem[0]);
  return qs_real(a0, a1, a2, a3);
}

static qd_real qd_from_qs(const qs_real &a) {
  return qd_real(static_cast<double>(a[0]), static_cast<double>(a[1]),
                 static_cast<double>(a[2]), static_cast<double>(a[3]));
}

/*======== Constants ========*/
namespace {
  qs_real make_qs(double a, double b, double c, double d = 0.0) {
    return qs_from_qd(qd_real(a, b, c, d));
  }
}

const qs_real qs_real::_2pi  = make_qs(6.283185307179586232e+00,  2.449293598294706414e-16, -5.989539619436679332e-33,  2.224908441726145253e-49);
const qs_real qs_real::_pi   = make_qs(3.141592653589793116e+00,  1.224646799147353207e-16, -2.994769809718339666e-33,  1.112454220863072626e-49);
const qs_real qs_real::_pi2  = make_qs(1.570796326794896558e+00,  6.123233995736766036e-17, -1.497384904859169833e-33,  5.562271104315363132e-50);
const qs_real qs_real::_pi4  = make_qs(7.853981633974482790e-01,  3.061616997868383018e-17, -7.486924524295849165e-34,  2.781135552157681566e-50);
const qs_real qs_real::_3pi4 = make_qs(2.356194490192344837e+00,  9.184850993605148e-17,     3.916898464750400e-33, -2.586798811020051e-49);
const qs_real qs_real::_e    = make_qs(2.718281828459045091e+00,  1.445646891729250158e-16, -2.127717108038176765e-33,  1.515630159841218973e-49);
const qs_real qs_real::_log2 = make_qs(6.931471805599452862e-01,  2.319046813846299558e-17,  5.707708438416212e-34,    -3.582432210601811e-50);
const qs_real qs_real::_log10= make_qs(2.302585092994045901e+00, -2.170756223382249351e-16, -9.984262454465776e-33, -4.023357454450206379e-49);

const qs_real qs_real::_nan = qs_real(qs::_f_nan, qs::_f_nan, qs::_f_nan, qs::_f_nan);
const qs_real qs_real::_inf = qs_real(qs::_f_inf, qs::_f_inf, qs::_f_inf, qs::_f_inf);

/* eps ~ 2^-93 ~ 1.01e-28 (4*24 - 3 working bits) */
const float qs_real::_eps = 1.0097420e-28f;
/* 2^(emin + 3p) = 2^(-126 + 72) = 2^-54 ~ 5.55112e-17 */
const float qs_real::_min_normalized = 5.55111512e-17f;
const qs_real qs_real::_max      = qs_real(3.4028234e+38f, 1.0141204e+31f, 3.0223145e+23f, 9.0094863e+15f);
const qs_real qs_real::_safe_max = qs_real(3.4028232e+38f, 1.0141204e+31f, 3.0223145e+23f, 9.0094863e+15f);
const int     qs_real::_ndigits  = 28;

/*======== Errors ========*/
void qs_real::error(const char *msg) {
  if (msg) cerr << "ERROR (qs_real): " << msg << endl;
}

/*======== sqrt ========*/
qs_real sqrt(const qs_real &a) {
  if (a.is_zero()) return qs_real(0.0f);
  if (a.is_negative()) {
    qs_real::error("sqrt: Negative argument.");
    return qs_real::_nan;
  }
  /* Newton iteration  r' = r + (0.5 - h*r^2)*r  for 1/sqrt(a), h = a/2.
     The first refinement runs in ds_real, the remaining two in qs_real,
     each as a pair of branch-free qw_fma calls.  This replaces the 0.0.2
     implementation, which promoted the argument to qd_real and ran three
     quad-double iterations. */
  ds_real ad(a[0], a[1]);
  ds_real hd = mul_pwr2(ad, 0.5f);
  ds_real rd = 1.0f / std::sqrt(a[0]);
  rd = dw_fma(dw_fma(-hd, sqr(rd), ds_real(0.5f)), rd, rd);

  ts_real at(a[0], a[1], a[2]);
  ts_real ht = mul_pwr2(at, 0.5f);
  ts_real rt(rd);
  rt = tw_fma(tw_fma(-ht, sqr(rt), ts_real(0.5f)), rt, rt);

  qs_real h = mul_pwr2(a, 0.5f);
  qs_real r(rt);
  r = qw_fma(qw_fma(-h, sqr(r), qs_real(0.5f)), r, r);

  r *= a;
  return r;
}

/* Reference (0.0.2) square root, kept for benchmarking. */
qs_real sqrt_legacy(const qs_real &a) {
  if (a.is_zero()) return qs_real(0.0f);
  if (a.is_negative()) {
    qs_real::error("sqrt_legacy: Negative argument.");
    return qs_real::_nan;
  }
  return qs_from_qd(sqrt(qd_from_qs(a)));
}

/*======== Powers ========*/
qs_real npwr(const qs_real &a, int n) {
  if (n == 0) {
    if (a.is_zero()) {
      qs_real::error("npwr: Invalid 0^0.");
      return qs_real::_nan;
    }
    return qs_real(1.0f);
  }
  qs_real r = a;
  qs_real s = 1.0f;
  int N = std::abs(n);
  if (N > 1) {
    while (N > 0) {
      if (N % 2 == 1) s *= r;
      N /= 2;
      if (N > 0) r = sqr(r);
    }
  } else {
    s = r;
  }
  if (n < 0) return qs_real(1.0f) / s;
  return s;
}

qs_real pow(const qs_real &a, int n)            { return npwr(a, n); }
qs_real pow(const qs_real &a, const qs_real &b) { return exp(b * log(a)); }

qs_real nroot(const qs_real &a, int n) {
  if (n <= 0) {
    qs_real::error("nroot: N must be positive.");
    return qs_real::_nan;
  }
  if (n % 2 == 0 && a.is_negative()) {
    qs_real::error("nroot: Negative argument.");
    return qs_real::_nan;
  }
  if (n == 1) return a;
  if (n == 2) return sqrt(a);
  if (a.is_zero()) return qs_real(0.0f);

  qs_real r = abs(a);
  qs_real x = std::exp(-std::log(static_cast<double>(r[0])) / n);
  for (int it = 0; it < 3; ++it) {
    x += x * (qs_real(1.0f) - r * npwr(x, n)) / static_cast<float>(n);
  }
  if (a[0] < 0.0f) x = -x;
  return qs_real(1.0f) / x;
}

/*======== Transcendentals (delegate to qd_real) ========*/
qs_real exp(const qs_real &a) { return qs_from_qd(exp(qd_from_qs(a))); }
qs_real log(const qs_real &a) {
  if (a.is_zero())     return -qs_real::_inf;
  if (a.is_negative()) { qs_real::error("log: Negative argument."); return qs_real::_nan; }
  return qs_from_qd(log(qd_from_qs(a)));
}
qs_real log10(const qs_real &a) { return log(a) / qs_real::_log10; }

qs_real sin(const qs_real &a) { return qs_from_qd(sin(qd_from_qs(a))); }
qs_real cos(const qs_real &a) { return qs_from_qd(cos(qd_from_qs(a))); }
qs_real tan(const qs_real &a) { return sin(a) / cos(a); }
void    sincos(const qs_real &a, qs_real &s, qs_real &c) {
  qd_real qs_, qc_;
  sincos(qd_from_qs(a), qs_, qc_);
  s = qs_from_qd(qs_);
  c = qs_from_qd(qc_);
}

/*======== Polyeval / Polyroot ========*/
qs_real polyeval(const qs_real *c, int n, const qs_real &x) {
  qs_real r = c[n];
  for (int i = n - 1; i >= 0; i--) r = r * x + c[i];
  return r;
}

qs_real polyroot(const qs_real *c, int n,
                 const qs_real &x0, int max_iter, float thresh) {
  qs_real x = x0;
  qs_real f, fp;
  qs_real *d = new qs_real[n];
  bool conv = false;
  float max_c = std::fabs(to_double(c[0]));

  if (thresh == 0.0f) thresh = qs_real::_eps;

  for (int i = 1; i <= n; i++) {
    float v = std::fabs(to_double(c[i]));
    if (v > max_c) max_c = v;
    d[i - 1] = c[i] * static_cast<float>(i);
  }
  thresh *= max_c;

  for (int i = 0; i < max_iter; i++) {
    f = polyeval(c, n, x);
    if (std::fabs(to_double(f)) < thresh) { conv = true; break; }
    fp = polyeval(d, n - 1, x);
    x -= (f / fp);
  }
  delete [] d;

  if (!conv) {
    qs_real::error("polyroot: Failed to converge.");
    return qs_real::_nan;
  }
  return x;
}

/*======== Random number ========*/
qs_real qsrand(void) {
  qs_real r = 0.0f;
  qs_real m = 1.0f;
  for (int i = 0; i < 8; i++) {
    m *= 3.0517578125e-05f; /* 2^-15 */
    float d = static_cast<float>(std::rand() % 32768);
    r += m * d;
  }
  return r;
}

/*======== String I/O ========*/
qs_real::qs_real(const char *s) {
  qd_real t;
  std::istringstream iss(s);
  iss >> t;
  *this = qs_from_qd(t);
}

qs_real &qs_real::operator=(const char *s) { *this = qs_real(s); return *this; }

string qs_real::to_string(int precision) const {
  return qd_from_qs(*this).to_string(precision);
}

ostream &operator<<(ostream &os, const qs_real &a) {
  std::ios_base::fmtflags fmt = os.flags();
  std::streamsize prec = os.precision();
  if (prec <= 0) prec = qs_real::_ndigits;
  return os << qd_from_qs(a).to_string(static_cast<int>(prec), 0, fmt);
}

istream &operator>>(istream &is, qs_real &a) {
  qd_real t;
  is >> t;
  a = qs_from_qd(t);
  return is;
}

/*======== ^ operator ========*/
qs_real qs_real::operator^(int n) const { return npwr(*this, n); }
