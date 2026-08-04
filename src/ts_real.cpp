/*
 * src/ts_real.cpp
 *
 * Implementation of the ts_real (triple-single, 3 x float) class.
 * Constants and transcendentals delegate to the higher-precision td_real
 * routines and the result is rounded back to ~72-bit (ts_real) precision.
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
#include <qd/ts_real.h>
#include <qd/td_real.h>

using std::cout;
using std::cerr;
using std::endl;
using std::ostream;
using std::istream;
using std::string;

/*======== Helpers: ts_real <-> td_real ========*/

/* td -> ts: extract the top ~72 bits as three floats.
   Use td_real arithmetic for the residual so that no bits are silently
   dropped between extractions.                                          */
static ts_real ts_from_td(const td_real &orig) {
  td_real rem = orig;
  float a0 = static_cast<float>(rem[0]);
  rem -= static_cast<double>(a0);
  float a1 = static_cast<float>(rem[0]);
  rem -= static_cast<double>(a1);
  float a2 = static_cast<float>(rem[0]);
  return ts_real(a0, a1, a2);
}

/* ts -> td: float values fit exactly in double. */
static td_real td_from_ts(const ts_real &a) {
  return td_real(static_cast<double>(a[0]),
                 static_cast<double>(a[1]),
                 static_cast<double>(a[2]));
}

/*======== Constants ========*/
namespace {
  ts_real make_ts(double hi, double mid, double lo = 0.0) {
    return ts_from_td(td_real(hi, mid, lo));
  }
}

const ts_real ts_real::_2pi  = make_ts(6.283185307179586232e+00,  2.449293598294706414e-16, -5.989539619436679332e-33);
const ts_real ts_real::_pi   = make_ts(3.141592653589793116e+00,  1.224646799147353207e-16, -2.994769809718339666e-33);
const ts_real ts_real::_pi2  = make_ts(1.570796326794896558e+00,  6.123233995736766036e-17, -1.497384904859169833e-33);
const ts_real ts_real::_pi4  = make_ts(7.853981633974482790e-01,  3.061616997868383018e-17, -7.486924524295849165e-34);
const ts_real ts_real::_3pi4 = make_ts(2.356194490192344837e+00,  9.184850993605148e-17,     3.916898464750400e-33);
const ts_real ts_real::_e    = make_ts(2.718281828459045091e+00,  1.445646891729250158e-16, -2.127717108038176765e-33);
const ts_real ts_real::_log2 = make_ts(6.931471805599452862e-01,  2.319046813846299558e-17,  5.707708438416212e-34);
const ts_real ts_real::_log10= make_ts(2.302585092994045901e+00, -2.170756223382249351e-16, -9.984262454465776e-33);

const ts_real ts_real::_nan = ts_real(qs::_f_nan, qs::_f_nan, qs::_f_nan);
const ts_real ts_real::_inf = ts_real(qs::_f_inf, qs::_f_inf, qs::_f_inf);

/* eps ~ 2^-69 ~ 1.69e-21 (3*24 - 3 working bits) */
const float ts_real::_eps = 1.6940659e-21f;
/* min normalized so that the smallest component stays normal:
   2^(emin + 2p) = 2^(-126 + 48) = 2^-78 ~ 3.30872e-24 */
const float ts_real::_min_normalized = 3.30872245e-24f;
const ts_real ts_real::_max      = ts_real(3.4028234e+38f, 1.0141204e+31f, 3.0223145e+23f);
const ts_real ts_real::_safe_max = ts_real(3.4028232e+38f, 1.0141204e+31f, 3.0223145e+23f);
const int     ts_real::_ndigits  = 21;

/*======== Errors ========*/
void ts_real::error(const char *msg) {
  if (msg) cerr << "ERROR (ts_real): " << msg << endl;
}

/*======== sqrt ========*/
ts_real sqrt(const ts_real &a) {
  if (a.is_zero()) return 0.0f;
  if (a.is_negative()) {
    ts_real::error("sqrt: Negative argument.");
    return ts_real::_nan;
  }
  /* Newton iteration  r' = r + (0.5 - h*r^2)*r  for 1/sqrt(a), h = a/2.
     The first refinement runs in ds_real (~14 digits is already plenty
     for the second one to reach the full ~21 digits of ts_real), and each
     step is a pair of branch-free tw_fma calls.  This replaces the 0.0.2
     implementation, which promoted the argument to td_real and ran three
     triple-double iterations. */
  ds_real ad(a[0], a[1]);
  ds_real hd = mul_pwr2(ad, 0.5f);
  ds_real rd = 1.0f / std::sqrt(a[0]);
  rd = dw_fma(dw_fma(-hd, sqr(rd), ds_real(0.5f)), rd, rd);

  ts_real h = mul_pwr2(a, 0.5f);
  ts_real r(rd);

  r = tw_fma(tw_fma(-h, sqr(r), ts_real(0.5f)), r, r);

  r *= a;
  return r;
}

/* Reference (0.0.2) square root, kept for benchmarking. */
ts_real sqrt_legacy(const ts_real &a) {
  if (a.is_zero()) return 0.0f;
  if (a.is_negative()) {
    ts_real::error("sqrt_legacy: Negative argument.");
    return ts_real::_nan;
  }
  return ts_from_td(sqrt(td_from_ts(a)));
}

/*======== Powers ========*/
ts_real npwr(const ts_real &a, int n) {
  if (n == 0) {
    if (a.is_zero()) {
      ts_real::error("npwr: Invalid 0^0.");
      return ts_real::_nan;
    }
    return ts_real(1.0f);
  }
  ts_real r = a;
  ts_real s = 1.0f;
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
  if (n < 0) return ts_real(1.0f) / s;
  return s;
}

ts_real pow(const ts_real &a, int n)            { return npwr(a, n); }
ts_real pow(const ts_real &a, const ts_real &b) { return exp(b * log(a)); }

ts_real nroot(const ts_real &a, int n) {
  if (n <= 0) {
    ts_real::error("nroot: N must be positive.");
    return ts_real::_nan;
  }
  if (n % 2 == 0 && a.is_negative()) {
    ts_real::error("nroot: Negative argument.");
    return ts_real::_nan;
  }
  if (n == 1) return a;
  if (n == 2) return sqrt(a);
  if (a.is_zero()) return ts_real(0.0f);

  ts_real r = abs(a);
  ts_real x = std::exp(-std::log(static_cast<double>(r[0])) / n);
  /* Newton iterations on f(x) = x^{-n} - a; do two passes for safety. */
  for (int it = 0; it < 2; ++it) {
    x += x * (1.0f - r * npwr(x, n)) / static_cast<float>(n);
  }
  if (a[0] < 0.0f) x = -x;
  return ts_real(1.0f) / x;
}

/*======== Transcendentals (delegate to td_real) ========*/
ts_real exp(const ts_real &a) { return ts_from_td(exp(td_from_ts(a))); }
ts_real log(const ts_real &a) {
  if (a.is_zero())     return -ts_real::_inf;
  if (a.is_negative()) { ts_real::error("log: Negative argument."); return ts_real::_nan; }
  return ts_from_td(log(td_from_ts(a)));
}
ts_real log10(const ts_real &a) { return log(a) / ts_real::_log10; }

ts_real sin(const ts_real &a) { return ts_from_td(sin(td_from_ts(a))); }
ts_real cos(const ts_real &a) { return ts_from_td(cos(td_from_ts(a))); }
ts_real tan(const ts_real &a) { return sin(a) / cos(a); }
void    sincos(const ts_real &a, ts_real &s, ts_real &c) {
  td_real ts, tc;
  sincos(td_from_ts(a), ts, tc);
  s = ts_from_td(ts);
  c = ts_from_td(tc);
}

/*======== Polyeval / Polyroot ========*/
ts_real polyeval(const ts_real *c, int n, const ts_real &x) {
  ts_real r = c[n];
  for (int i = n - 1; i >= 0; i--) r = r * x + c[i];
  return r;
}

ts_real polyroot(const ts_real *c, int n,
                 const ts_real &x0, int max_iter, float thresh) {
  ts_real x = x0;
  ts_real f, fp;
  ts_real *d = new ts_real[n];
  bool conv = false;
  float max_c = std::fabs(to_double(c[0]));

  if (thresh == 0.0f) thresh = ts_real::_eps;

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
    ts_real::error("polyroot: Failed to converge.");
    return ts_real::_nan;
  }
  return x;
}

/*======== Random number ========*/
ts_real tsrand(void) {
  ts_real r = 0.0f;
  ts_real m = 1.0f;
  for (int i = 0; i < 6; i++) {
    m *= 3.0517578125e-05f; /* 2^-15 */
    float d = static_cast<float>(std::rand() % 32768);
    r += m * d;
  }
  return r;
}

/*======== String I/O ========*/
ts_real::ts_real(const char *s) {
  td_real t;
  std::istringstream iss(s);
  iss >> t;
  *this = ts_from_td(t);
}

ts_real &ts_real::operator=(const char *s) { *this = ts_real(s); return *this; }

string ts_real::to_string(int precision) const {
  return td_from_ts(*this).to_string(precision);
}

ostream &operator<<(ostream &os, const ts_real &a) {
  std::ios_base::fmtflags fmt = os.flags();
  std::streamsize prec = os.precision();
  if (prec <= 0) prec = ts_real::_ndigits;
  return os << td_from_ts(a).to_string(static_cast<int>(prec), 0, fmt);
}

istream &operator>>(istream &is, ts_real &a) {
  td_real t;
  is >> t;
  a = ts_from_td(t);
  return is;
}

/*======== ^ operator ========*/
ts_real ts_real::operator^(int n) const { return npwr(*this, n); }
