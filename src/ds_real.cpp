/*
 * src/ds_real.cpp
 *
 * Implementation of the ds_real (single-double, 2 x float) class.
 * Constants and transcendentals delegate to the higher-precision dd_real
 * routines and the result is rounded back to ~48-bit (ds_real) precision.
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
#include <qd/ds_real.h>
#include <qd/dd_real.h>

using std::cout;
using std::cerr;
using std::endl;
using std::ostream;
using std::istream;
using std::string;

/*======== Helpers: convert between ds_real and dd_real ========*/

/* Convert a dd_real value back into a ds_real, keeping ~48 bits of mantissa.
   Use dd_real arithmetic for the residual so we don't lose any bits.  */
static ds_real ds_from_dd(const dd_real &orig) {
  dd_real rem = orig;
  float a0 = static_cast<float>(rem._hi());
  rem -= static_cast<double>(a0);
  float a1 = static_cast<float>(rem._hi());
  return ds_real(a0, a1);
}

static dd_real dd_from_ds(const ds_real &a) {
  /* float values fit exactly in double; just place them. */
  return dd_real(static_cast<double>(a.x[0]), static_cast<double>(a.x[1]));
}

/*======== Constants ========*/

namespace {
  /* Build a ds_real from a (hi,lo) double-double approximation. */
  ds_real make_ds(double hi, double lo) {
    return ds_from_dd(dd_real(hi, lo));
  }
}

const ds_real ds_real::_2pi  = make_ds(6.283185307179586232e+00,  2.449293598294706414e-16);
const ds_real ds_real::_pi   = make_ds(3.141592653589793116e+00,  1.224646799147353207e-16);
const ds_real ds_real::_pi2  = make_ds(1.570796326794896558e+00,  6.123233995736766036e-17);
const ds_real ds_real::_pi4  = make_ds(7.853981633974482790e-01,  3.061616997868383018e-17);
const ds_real ds_real::_3pi4 = make_ds(2.356194490192344837e+00,  9.184850993605148e-17);
const ds_real ds_real::_e    = make_ds(2.718281828459045091e+00,  1.445646891729250158e-16);
const ds_real ds_real::_log2 = make_ds(6.931471805599452862e-01,  2.319046813846299558e-17);
const ds_real ds_real::_log10= make_ds(2.302585092994045901e+00, -2.170756223382249351e-16);

const ds_real ds_real::_nan = ds_real(qs::_f_nan, qs::_f_nan);
const ds_real ds_real::_inf = ds_real(qs::_f_inf, qs::_f_inf);

/* eps = 2^-46 ~ 1.42e-14    (2*24 - 2 working bits) */
const float ds_real::_eps = 1.4210854715202004e-14f;
/* float min normalized = 2^-126 ~ 1.175e-38;
   for ds we need room for the secondary float -> use 2^-126 + 24 bits */
const float ds_real::_min_normalized = 1.9721523e-31f;        /* ~ 2^-102 */
const ds_real ds_real::_max       = ds_real(3.4028234e+38f, 1.0141204e+31f);
const ds_real ds_real::_safe_max  = ds_real(3.4028232e+38f, 1.0141204e+31f);
const int     ds_real::_ndigits   = 14;

/*======== Error / errors ========*/
void ds_real::error(const char *msg) {
  if (msg) cerr << "ERROR (ds_real): " << msg << endl;
}

/*======== Float division -> ds_real ========*/
/* Defined inline in header. */

/*======== sqrt ========*/
ds_real sqrt(const ds_real &a) {
  if (a.is_zero()) return 0.0f;
  if (a.is_negative()) {
    ds_real::error("sqrt: Negative argument.");
    return ds_real::_nan;
  }
  /* Karp's trick. */
  float x = 1.0f / std::sqrt(a.x[0]);
  float ax = a.x[0] * x;
  return ds_real::add(ax, (a - ds_real::sqr(ax)).x[0] * (x * 0.5f));
}

ds_real ds_real::sqrt(float a) {
  return ::sqrt(ds_real(a));
}

/*======== Powers ========*/
ds_real npwr(const ds_real &a, int n) {
  if (n == 0) {
    if (a.is_zero()) {
      ds_real::error("npwr: Invalid 0^0.");
      return ds_real::_nan;
    }
    return 1.0f;
  }
  ds_real r = a;
  ds_real s = 1.0f;
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
  if (n < 0) return 1.0f / s;
  return s;
}

ds_real pow(const ds_real &a, int n)             { return npwr(a, n); }
ds_real pow(const ds_real &a, const ds_real &b)  { return exp(b * log(a)); }

ds_real nroot(const ds_real &a, int n) {
  if (n <= 0) {
    ds_real::error("nroot: N must be positive.");
    return ds_real::_nan;
  }
  if (n % 2 == 0 && a.is_negative()) {
    ds_real::error("nroot: Negative argument.");
    return ds_real::_nan;
  }
  if (n == 1) return a;
  if (n == 2) return sqrt(a);
  if (a.is_zero()) return 0.0f;

  ds_real r = abs(a);
  ds_real x = std::exp(-std::log(static_cast<double>(r.x[0])) / n);
  /* Newton's iteration on f(x) = x^{-n} - a */
  x += x * (1.0f - r * npwr(x, n)) / static_cast<float>(n);
  if (a.x[0] < 0.0f) x = -x;
  return 1.0f / x;
}

/*======== Transcendentals (delegate to dd_real) ========*/
ds_real exp(const ds_real &a)  { return ds_from_dd(exp (dd_from_ds(a))); }
ds_real log(const ds_real &a)  {
  if (a.is_zero())     return -ds_real::_inf;
  if (a.is_negative()) { ds_real::error("log: Negative argument."); return ds_real::_nan; }
  return ds_from_dd(log(dd_from_ds(a)));
}
ds_real log10(const ds_real &a){ return log(a) / ds_real::_log10; }

ds_real sin(const ds_real &a) { return ds_from_dd(sin(dd_from_ds(a))); }
ds_real cos(const ds_real &a) { return ds_from_dd(cos(dd_from_ds(a))); }
ds_real tan(const ds_real &a) { return sin(a) / cos(a); }
void    sincos(const ds_real &a, ds_real &s, ds_real &c) {
  dd_real ds, dc;
  sincos(dd_from_ds(a), ds, dc);
  s = ds_from_dd(ds);
  c = ds_from_dd(dc);
}

/*======== Polyeval / Polyroot ========*/
ds_real polyeval(const ds_real *c, int n, const ds_real &x) {
  ds_real r = c[n];
  for (int i = n - 1; i >= 0; i--) {
    r = r * x + c[i];
  }
  return r;
}

ds_real polyroot(const ds_real *c, int n,
                 const ds_real &x0, int max_iter, float thresh) {
  ds_real x = x0;
  ds_real f, fp;
  ds_real *d = new ds_real[n];
  bool conv = false;
  float max_c = std::fabs(to_double(c[0]));

  if (thresh == 0.0f) thresh = ds_real::_eps;

  /* Compute coefficients of derivative */
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
    ds_real::error("polyroot: Failed to converge.");
    return ds_real::_nan;
  }
  return x;
}

/*======== Random number ========*/
ds_real dsrand(void) {
  ds_real r = 0.0f;
  ds_real m = 1.0f;
  for (int i = 0; i < 4; i++) {
    m *= 3.0517578125e-05f; /* 2^-15 */
    float d = static_cast<float>(std::rand() % 32768);
    r += m * d;
  }
  return r;
}

/*======== String parsing/printing ========*/
ds_real::ds_real(const char *s) {
  dd_real t;
  std::istringstream iss(s);
  iss >> t;
  *this = ds_from_dd(t);
}

ds_real &ds_real::operator=(const char *s) { *this = ds_real(s); return *this; }

string ds_real::to_string(int precision) const {
  return dd_from_ds(*this).to_string(precision);
}

ostream &operator<<(ostream &os, const ds_real &a) {
  std::ios_base::fmtflags fmt = os.flags();
  std::streamsize prec = os.precision();
  if (prec <= 0) prec = ds_real::_ndigits;
  /* Delegate stream formatting to dd_real for richer behavior. */
  return os << dd_from_ds(a).to_string(static_cast<int>(prec), 0, fmt);
}

istream &operator>>(istream &is, ds_real &a) {
  dd_real t;
  is >> t;
  a = ds_from_dd(t);
  return is;
}

/*======== ^ operator ========*/
ds_real ds_real::operator^(int n) const { return npwr(*this, n); }
