/*
 * include/qd/ts_real.h
 *
 * Triple-single precision (>= 72-bit significand) floating point arithmetic
 * package.  A ts_real is represented by an unevaluated sum of three IEEE
 * single-precision floats x[0] + x[1] + x[2] in non-overlapping form,
 * giving roughly 21 decimal digits of precision.
 *
 * Modeled directly on td_real / qd_real (LBNL).
 *
 * Created for the dtq package extension (2026).
 */
#ifndef _QD_TS_REAL_H
#define _QD_TS_REAL_H

#include <iostream>
#include <string>
#include <limits>
#include <qd/qd_config.h>
#include <qd/fs_inline.h>
#include <qd/ds_real.h>

/* Full ~72-bit precision holds only while the lowest limb (~|x|*2^-48)
   stays in the float normal range, i.e. |x| >= 2^-78; below that, precision
   degrades by about one bit per halving (representation limit of the
   format itself, independent of the algorithm). */
struct QD_API ts_real {
  float x[3];

  void renorm();
  void renorm(float &e);

  ts_real(float x0, float x1, float x2) { x[0] = x0; x[1] = x1; x[2] = x2; }
  explicit ts_real(const float *xx) { x[0] = xx[0]; x[1] = xx[1]; x[2] = xx[2]; }

  static const ts_real _2pi;
  static const ts_real _pi;
  static const ts_real _3pi4;
  static const ts_real _pi2;
  static const ts_real _pi4;
  static const ts_real _e;
  static const ts_real _log2;
  static const ts_real _log10;
  static const ts_real _nan;
  static const ts_real _inf;

  static const float _eps;
  static const float _min_normalized;
  static const ts_real _max;
  static const ts_real _safe_max;
  static const int _ndigits;

  ts_real() { x[0] = x[1] = x[2] = 0.0f; }
  ts_real(const char *s);
  ts_real(const ds_real &a) { x[0] = a._hi(); x[1] = a._lo(); x[2] = 0.0f; }
  ts_real(float d)  { x[0] = d; x[1] = x[2] = 0.0f; }
  ts_real(double d) { x[0] = static_cast<float>(d); x[1] = x[2] = 0.0f; }
  ts_real(int i)    { x[0] = static_cast<float>(i); x[1] = x[2] = 0.0f; }

  float operator[](int i) const { return x[i]; }
  float &operator[](int i) { return x[i]; }

  static void error(const char *msg);

  bool isnan() const { return QD_ISNAN(x[0]) || QD_ISNAN(x[1]) || QD_ISNAN(x[2]); }
  bool isfinite() const { return QD_ISFINITE(x[0]); }
  bool isinf() const { return QD_ISINF(x[0]); }

  static ts_real ieee_add(const ts_real &a, const ts_real &b);
  static ts_real sloppy_add(const ts_real &a, const ts_real &b);

  ts_real &operator+=(float a);
  ts_real &operator+=(const ds_real &a);
  ts_real &operator+=(const ts_real &a);

  ts_real &operator-=(float a);
  ts_real &operator-=(const ds_real &a);
  ts_real &operator-=(const ts_real &a);

  static ts_real sloppy_mul(const ts_real &a, const ts_real &b);
  static ts_real accurate_mul(const ts_real &a, const ts_real &b);

  ts_real &operator*=(float a);
  ts_real &operator*=(const ds_real &a);
  ts_real &operator*=(const ts_real &a);

  static ts_real sloppy_div(const ts_real &a, const ts_real &b);
  static ts_real accurate_div(const ts_real &a, const ts_real &b);
  /* Division built on the branch-free triple-word FMA (tw_fma). */
  static ts_real fma_div(const ts_real &a, const ts_real &b);

  ts_real &operator/=(float a);
  ts_real &operator/=(const ts_real &a);

  ts_real operator^(int n) const;
  ts_real operator-() const { return ts_real(-x[0], -x[1], -x[2]); }

  ts_real &operator=(float a)  { x[0] = a; x[1] = x[2] = 0.0f; return *this; }
  ts_real &operator=(double a) { x[0] = static_cast<float>(a); x[1] = x[2] = 0.0f; return *this; }
  ts_real &operator=(int a)    { x[0] = static_cast<float>(a); x[1] = x[2] = 0.0f; return *this; }
  ts_real &operator=(const ds_real &a) { x[0] = a._hi(); x[1] = a._lo(); x[2] = 0.0f; return *this; }
  ts_real &operator=(const char *s);

  bool is_zero()     const { return x[0] == 0.0f; }
  bool is_one()      const { return x[0] == 1.0f && x[1] == 0.0f && x[2] == 0.0f; }
  bool is_positive() const { return x[0] >  0.0f; }
  bool is_negative() const { return x[0] <  0.0f; }

  static ts_real rand(void);

  std::string to_string(int precision = _ndigits) const;
};

namespace std {
  template <>
  class numeric_limits<ts_real> : public numeric_limits<float> {
  public:
    inline static float epsilon() { return ts_real::_eps; }
    inline static float min() { return ts_real::_min_normalized; }
    inline static ts_real max() { return ts_real::_max; }
    inline static ts_real safe_max() { return ts_real::_safe_max; }
    static const int digits = 72;
    static const int digits10 = 21;
  };
}

QD_API ts_real polyeval(const ts_real *c, int n, const ts_real &x);
QD_API ts_real polyroot(const ts_real *c, int n,
    const ts_real &x0, int max_iter = 48, float thresh = 0.0f);

QD_API ts_real tsrand(void);
QD_API ts_real sqrt(const ts_real &a);

QD_API inline bool isnan(const ts_real &a) { return a.isnan(); }
QD_API inline bool isfinite(const ts_real &a) { return a.isfinite(); }
QD_API inline bool isinf(const ts_real &a) { return a.isinf(); }

QD_API ts_real mul_pwr2(const ts_real &a, float b);

/*********** Branch-free triple-word FMA (T. Kouya) ************/
/* tw_fma(a, b, c) computes  a * b + c  as a single fused triple-word
   operation, without renormalizing a * b on its own. */
QD_API ts_real tw_fma(const ts_real &a, const ts_real &b, const ts_real &c);
QD_API ts_real tw_fma(const ts_real &a, float b, const ts_real &c);

/* Generic spelling; same operation. */
QD_API ts_real fma(const ts_real &a, const ts_real &b, const ts_real &c);
QD_API ts_real fma(const ts_real &a, float b, const ts_real &c);

/* Reference (pre-0.0.3) implementation, kept for benchmarking. */
QD_API ts_real sqrt_legacy(const ts_real &a);

QD_API ts_real operator+(const ts_real &a, const ts_real &b);
QD_API ts_real operator+(const ds_real &a, const ts_real &b);
QD_API ts_real operator+(const ts_real &a, const ds_real &b);
QD_API ts_real operator+(const ts_real &a, float b);
QD_API ts_real operator+(float a, const ts_real &b);

QD_API ts_real operator-(const ts_real &a, const ts_real &b);
QD_API ts_real operator-(const ts_real &a, float b);
QD_API ts_real operator-(float a, const ts_real &b);

QD_API ts_real operator*(const ts_real &a, const ts_real &b);
QD_API ts_real operator*(const ts_real &a, float b);
QD_API ts_real operator*(float a, const ts_real &b);

QD_API ts_real operator/(const ts_real &a, const ts_real &b);
QD_API ts_real operator/(const ts_real &a, float b);
QD_API ts_real operator/(float a, const ts_real &b);

QD_API ts_real sqr(const ts_real &a);
QD_API ts_real sqrt(const ts_real &a);
QD_API ts_real pow(const ts_real &a, int n);
QD_API ts_real pow(const ts_real &a, const ts_real &b);
QD_API ts_real npwr(const ts_real &a, int n);
QD_API ts_real nroot(const ts_real &a, int n);

ds_real to_ds_real(const ts_real &a);
QD_API float to_double(const ts_real &a);
QD_API float to_float (const ts_real &a);
QD_API int   to_int   (const ts_real &a);

QD_API bool operator==(const ts_real &a, const ts_real &b);
QD_API bool operator==(const ts_real &a, float b);
QD_API bool operator==(float a, const ts_real &b);

QD_API bool operator<(const ts_real &a, const ts_real &b);
QD_API bool operator<(const ts_real &a, float b);
QD_API bool operator<(float a, const ts_real &b);

QD_API bool operator>(const ts_real &a, const ts_real &b);
QD_API bool operator>(const ts_real &a, float b);
QD_API bool operator>(float a, const ts_real &b);

QD_API bool operator<=(const ts_real &a, const ts_real &b);
QD_API bool operator<=(const ts_real &a, float b);
QD_API bool operator<=(float a, const ts_real &b);

QD_API bool operator>=(const ts_real &a, const ts_real &b);
QD_API bool operator>=(const ts_real &a, float b);
QD_API bool operator>=(float a, const ts_real &b);

QD_API bool operator!=(const ts_real &a, const ts_real &b);
QD_API bool operator!=(const ts_real &a, float b);
QD_API bool operator!=(float a, const ts_real &b);

QD_API ts_real fabs(const ts_real &a);
QD_API ts_real abs(const ts_real &a);

QD_API ts_real ldexp(const ts_real &a, int n);
QD_API ts_real nint(const ts_real &a);
QD_API ts_real floor(const ts_real &a);
QD_API ts_real ceil(const ts_real &a);
QD_API ts_real aint(const ts_real &a);

QD_API ts_real exp(const ts_real &a);
QD_API ts_real log(const ts_real &a);
QD_API ts_real log10(const ts_real &a);

QD_API ts_real sin(const ts_real &a);
QD_API ts_real cos(const ts_real &a);
QD_API ts_real tan(const ts_real &a);
QD_API void    sincos(const ts_real &a, ts_real &s, ts_real &c);

QD_API std::ostream &operator<<(std::ostream &s, const ts_real &a);
QD_API std::istream &operator>>(std::istream &s, ts_real &a);

/********** Renormalization (in qs::) **********/
namespace qs {
inline void renorm3(float &c0, float &c1, float &c2) {
  float s0, s1, s2 = 0.0f;
  if (QD_ISINF(c0)) return;
  s0 = qs::quick_two_sum(c1, c2, c2);
  c0 = qs::quick_two_sum(c0, s0, c1);
  s0 = c0; s1 = c1;
  if (s1 != 0.0f) {
    s1 = qs::quick_two_sum(s1, c2, s2);
  } else {
    s0 = qs::quick_two_sum(s0, c2, s1);
  }
  c0 = s0; c1 = s1; c2 = s2;
}
inline void renorm3(float &c0, float &c1, float &c2, float &c3) {
  float s0, s1, s2 = 0.0f;
  if (QD_ISINF(c0)) return;
  s0 = qs::quick_two_sum(c2, c3, c3);
  s0 = qs::quick_two_sum(c1, s0, c2);
  c0 = qs::quick_two_sum(c0, s0, c1);
  s0 = c0; s1 = c1;
  if (s1 != 0.0f) {
    s1 = qs::quick_two_sum(s1, c2, s2);
    if (s2 != 0.0f) s2 += c3;
    else            s1 = qs::quick_two_sum(s1, c3, s2);
  } else {
    s0 = qs::quick_two_sum(s0, c2, s1);
    if (s1 != 0.0f) s1 = qs::quick_two_sum(s1, c3, s2);
    else            s0 = qs::quick_two_sum(s0, c3, s1);
  }
  c0 = s0; c1 = s1; c2 = s2;
}
}

inline void ts_real::renorm()             { qs::renorm3(x[0], x[1], x[2]); }
inline void ts_real::renorm(float &e)     { qs::renorm3(x[0], x[1], x[2], e); }

/********** Inline arithmetic **********/

inline ts_real operator+(const ts_real &a, float b) {
  float c0, c1, c2, e;
  c0 = qs::two_sum(a[0], b, e);
  c1 = qs::two_sum(a[1], e, e);
  c2 = qs::two_sum(a[2], e, e);
  qs::renorm3(c0, c1, c2, e);
  return ts_real(c0, c1, c2);
}
inline ts_real operator+(float a, const ts_real &b) { return b + a; }

inline ts_real ts_real::sloppy_add(const ts_real &a, const ts_real &b) {
  float s0, s1, s2;
  float t0, t1, t2;
  s0 = qs::two_sum(a[0], b[0], t0);
  s1 = qs::two_sum(a[1], b[1], t1);
  s2 = qs::two_sum(a[2], b[2], t2);
  s1 = qs::two_sum(s1, t0, t0);
  qs::three_sum(s2, t0, t1);
  float tail = t0 + t1 + t2;
  qs::renorm3(s0, s1, s2, tail);
  return ts_real(s0, s1, s2);
}

inline ts_real operator+(const ts_real &a, const ts_real &b) {
  return ts_real::sloppy_add(a, b);
}
inline ts_real operator+(const ts_real &a, const ds_real &b) {
  return a + ts_real(b);
}
inline ts_real operator+(const ds_real &a, const ts_real &b) {
  return ts_real(a) + b;
}

inline ts_real operator-(const ts_real &a, const ts_real &b) { return a + (-b); }
inline ts_real operator-(const ts_real &a, float b)          { return a + (-b); }
inline ts_real operator-(float a, const ts_real &b)          { return a + (-b); }

inline ts_real operator*(const ts_real &a, float b) {
  float p0, p1, p2;
  float q0, q1;
  p0 = qs::two_prod(a[0], b, q0);
  p1 = qs::two_prod(a[1], b, q1);
  p2 = a[2] * b;
  qs::three_sum(p1, p2, q0);
  p2 = p2 + q0 + q1;
  qs::renorm3(p0, p1, p2);
  return ts_real(p0, p1, p2);
}
inline ts_real operator*(float a, const ts_real &b) { return b * a; }

inline ts_real ts_real::sloppy_mul(const ts_real &a, const ts_real &b) {
  float p0, p1, p2;
  float q0, q1, q2;
  p0 = qs::two_prod(a[0], b[0], q0);
  p1 = qs::two_prod(a[0], b[1], q1);
  p2 = qs::two_prod(a[1], b[0], q2);
  qs::three_sum(p1, p2, q0);
  p2 = p2 + q1 + q2;
  /* lower-order cross terms */
  p2 += a[0] * b[2] + a[2] * b[0] + a[1] * b[1];
  qs::renorm3(p0, p1, p2);
  return ts_real(p0, p1, p2);
}
inline ts_real operator*(const ts_real &a, const ts_real &b) {
  return ts_real::sloppy_mul(a, b);
}

inline ts_real mul_pwr2(const ts_real &a, float b) {
  return ts_real(a[0]*b, a[1]*b, a[2]*b);
}
inline ts_real ldexp(const ts_real &a, int n) {
  return ts_real(std::ldexp(a[0], n), std::ldexp(a[1], n), std::ldexp(a[2], n));
}

// 2026-08-04 T.Kouya
// Branch free algorithm: triple-word FMA,  z = a * b + c  (float limbs).
// Same network as the td_real version; see include/qd/td_inline.h.
/* TW-FMA  z = a * b + c   (72 flops)
   Machine-proved with FPANVerifier + z3 5.0.0 (ACS2026 formulation):
     error bound      |z-(ab+c)| <= 184 u^3 (|ab|+|c|)
     every FastTwoSum precondition  exp(x) >= exp(y)
     non-overlapping output         z0 |> z1 |> z2   (strongly_dominates)
   The precision p is symbolic, so the proofs hold for binary32/64/128 alike.
   Normalization repeats a cascade over adjacent pairs; the pass count is the
   smallest for which the non-overlap is provable (DW 1 / TW 3 / QW 5). */
inline ts_real tw_fma(const ts_real &a, const ts_real &b, const ts_real &c) {
  float P00, E00, P01, E01, P10, E10, P02, P11, P20, sg, G;
  float A, q1, q2, q3, B, r, m1, m2, w0, w1, w2;
  P00 = qs::two_prod(a[0], b[0], E00);
  P01 = qs::two_prod(a[0], b[1], E01);
  P10 = qs::two_prod(a[1], b[0], E10);
  P02 = a[0] * b[2];
  P11 = a[1] * b[1];
  P20 = a[2] * b[0];
  sg  = (P02 + P20) + P11;
  G   = ((E01 + E10) + sg) + c[2];
  A   = qs::two_sum(P01, P10, q1);
  A   = qs::two_sum(A, E00, q2);
  A   = qs::two_sum(A, c[1], q3);
  G   = G + ((q1 + q2) + q3);
  B   = qs::two_sum(P00, c[0], r);
  m1  = qs::two_sum(r, A, m2);
  m2  = m2 + G;
  float z0, z1, z2;
  /* renormalization pass 1/3 */
  w0 = qs::quick_two_sum(B, m1, w1);
  w1 = qs::two_sum(w1, m2, w2);
  /* renormalization pass 2/3 */
  w0 = qs::two_sum(w0, w1, w1);
  w1 = qs::quick_two_sum(w1, w2, w2);
  /* renormalization pass 3/3 */
  z0 = qs::quick_two_sum(w0, w1, w1);
  z1 = qs::quick_two_sum(w1, w2, z2);
  return ts_real(z0, z1, z2);
}

/* TW-FMA  z = a * b + c   (72 flops, scalar multiplier)
   Machine-proved with FPANVerifier + z3 5.0.0 (ACS2026 formulation):
     error bound      |z-(ab+c)| <= 184 u^3 (|ab|+|c|)
     every FastTwoSum precondition  exp(x) >= exp(y)
     non-overlapping output         z0 |> z1 |> z2   (strongly_dominates)
   The precision p is symbolic, so the proofs hold for binary32/64/128 alike.
   Normalization repeats a cascade over adjacent pairs; the pass count is the
   smallest for which the non-overlap is provable (DW 1 / TW 3 / QW 5). */

/* div/sqrt-safe variant (84 flops): the Newton iterations of division and
   square root receive residuals that are NOT non-overlapping expansions, so no
   FastTwoSum precondition can be claimed for any gate.  A FastTwoSum whose
   precondition fails does not even satisfy s+e=a+b -- it is no longer an EFT --
   so this variant uses TwoSum everywhere, with the same pass count as the
   standard variant.  A safe variant is never cheaper than the standard one. */
inline ts_real tw_fma_safe(const ts_real &a, float b, const ts_real &c) {
  float P00, E00, P01, E01, P10, E10, P02, P11, P20, sg, G;
  float A, q1, q2, q3, B, r, m1, m2, w0, w1, w2;
  P00 = qs::two_prod(a[0], b, E00);
  P01 = 0.0; E01 = 0.0;
  P10 = qs::two_prod(a[1], b, E10);
  P02 = 0.0;
  P11 = 0.0;
  P20 = a[2] * b;
  sg  = (P02 + P20) + P11;
  G   = ((E01 + E10) + sg) + c[2];
  A   = qs::two_sum(P01, P10, q1);
  A   = qs::two_sum(A, E00, q2);
  A   = qs::two_sum(A, c[1], q3);
  G   = G + ((q1 + q2) + q3);
  B   = qs::two_sum(P00, c[0], r);
  m1  = qs::two_sum(r, A, m2);
  m2  = m2 + G;
  float z0, z1, z2;
  /* renormalization pass 1/3 */
  w0 = qs::two_sum(B, m1, w1);
  w1 = qs::two_sum(w1, m2, w2);
  /* renormalization pass 2/3 */
  w0 = qs::two_sum(w0, w1, w1);
  w1 = qs::two_sum(w1, w2, w2);
  /* renormalization pass 3/3 */
  z0 = qs::two_sum(w0, w1, w1);
  z1 = qs::two_sum(w1, w2, z2);
  return ts_real(z0, z1, z2);
}

inline ts_real tw_fma(const ts_real &a, float b, const ts_real &c) {
  float P00, E00, P01, E01, P10, E10, P02, P11, P20, sg, G;
  float A, q1, q2, q3, B, r, m1, m2, w0, w1, w2;
  P00 = qs::two_prod(a[0], b, E00);
  P01 = 0.0; E01 = 0.0;
  P10 = qs::two_prod(a[1], b, E10);
  P02 = 0.0;
  P11 = 0.0;
  P20 = a[2] * b;
  sg  = (P02 + P20) + P11;
  G   = ((E01 + E10) + sg) + c[2];
  A   = qs::two_sum(P01, P10, q1);
  A   = qs::two_sum(A, E00, q2);
  A   = qs::two_sum(A, c[1], q3);
  G   = G + ((q1 + q2) + q3);
  B   = qs::two_sum(P00, c[0], r);
  m1  = qs::two_sum(r, A, m2);
  m2  = m2 + G;
  float z0, z1, z2;
  /* renormalization pass 1/3 */
  w0 = qs::quick_two_sum(B, m1, w1);
  w1 = qs::two_sum(w1, m2, w2);
  /* renormalization pass 2/3 */
  w0 = qs::two_sum(w0, w1, w1);
  w1 = qs::quick_two_sum(w1, w2, w2);
  /* renormalization pass 3/3 */
  z0 = qs::quick_two_sum(w0, w1, w1);
  z1 = qs::quick_two_sum(w1, w2, z2);
  return ts_real(z0, z1, z2);
}

inline ts_real fma(const ts_real &a, const ts_real &b, const ts_real &c) {
  return tw_fma(a, b, c);
}

inline ts_real fma(const ts_real &a, float b, const ts_real &c) {
  return tw_fma(a, b, c);
}

inline ts_real ts_real::sloppy_div(const ts_real &a, const ts_real &b) {
  float q0, q1, q2;
  ts_real r;
  q0 = a[0] / b[0];
  r = a - b * q0;
  q1 = r[0] / b[0];
  r = r - b * q1;
  q2 = r[0] / b[0];
  qs::renorm3(q0, q1, q2);
  return ts_real(q0, q1, q2);
}
inline ts_real ts_real::accurate_div(const ts_real &a, const ts_real &b) {
  float q0, q1, q2, q3;
  ts_real r;
  q0 = a[0] / b[0];
  r = a - b * q0;
  q1 = r[0] / b[0];
  r = r - b * q1;
  q2 = r[0] / b[0];
  r = r - b * q2;
  q3 = r[0] / b[0];
  qs::renorm3(q0, q1, q2, q3);
  return ts_real(q0, q1, q2);
}
/* Long division driven by the branch-free triple-word FMA.  Same
   correction sequence as sloppy_div (the ts_real default), but each
   residual  r <- r - q * b  is one fused tw_fma. */
inline ts_real ts_real::fma_div(const ts_real &a, const ts_real &b) {
  float q0, q1, q2;
  ts_real r;

  q0 = a[0] / b[0];
  r = tw_fma(b, -q0, a);
  q1 = r[0] / b[0];
  r = tw_fma(b, -q1, r);
  q2 = r[0] / b[0];

  qs::renorm3(q0, q1, q2);
  return ts_real(q0, q1, q2);
}

inline ts_real operator/(const ts_real &a, const ts_real &b) {
#ifdef QD_NO_FMA_DIV
  return ts_real::sloppy_div(a, b);
#else
  return ts_real::fma_div(a, b);
#endif
}
inline ts_real operator/(const ts_real &a, float b) { return a / ts_real(b); }
inline ts_real operator/(float a, const ts_real &b) { return ts_real(a) / b; }

inline ts_real &ts_real::operator+=(float a)         { *this = *this + a; return *this; }
inline ts_real &ts_real::operator+=(const ds_real &a){ *this = *this + a; return *this; }
inline ts_real &ts_real::operator+=(const ts_real &a){ *this = *this + a; return *this; }
inline ts_real &ts_real::operator-=(float a)         { *this = *this - a; return *this; }
inline ts_real &ts_real::operator-=(const ds_real &a){ *this = *this - ts_real(a); return *this; }
inline ts_real &ts_real::operator-=(const ts_real &a){ *this = *this - a; return *this; }
inline ts_real &ts_real::operator*=(float a)         { *this = *this * a; return *this; }
inline ts_real &ts_real::operator*=(const ds_real &a){ *this = *this * ts_real(a); return *this; }
inline ts_real &ts_real::operator*=(const ts_real &a){ *this = *this * a; return *this; }
inline ts_real &ts_real::operator/=(float a)         { *this = *this / a; return *this; }
inline ts_real &ts_real::operator/=(const ts_real &a){ *this = *this / a; return *this; }

/********** Comparisons **********/
inline bool operator==(const ts_real &a, const ts_real &b) {
  return a[0] == b[0] && a[1] == b[1] && a[2] == b[2];
}
inline bool operator==(const ts_real &a, float b) { return a[0] == b && a[1] == 0.0f && a[2] == 0.0f; }
inline bool operator==(float a, const ts_real &b) { return b == a; }
inline bool operator!=(const ts_real &a, const ts_real &b) { return !(a == b); }
inline bool operator!=(const ts_real &a, float b)          { return !(a == b); }
inline bool operator!=(float a, const ts_real &b)          { return !(a == b); }

inline bool operator<(const ts_real &a, const ts_real &b) {
  if (a[0] != b[0]) return a[0] < b[0];
  if (a[1] != b[1]) return a[1] < b[1];
  return a[2] < b[2];
}
inline bool operator<(const ts_real &a, float b) {
  if (a[0] != b) return a[0] < b;
  if (a[1] != 0.0f) return a[1] < 0.0f;
  return a[2] < 0.0f;
}
inline bool operator<(float a, const ts_real &b) {
  if (a != b[0]) return a < b[0];
  if (b[1] != 0.0f) return 0.0f < b[1];
  return 0.0f < b[2];
}

inline bool operator>(const ts_real &a, const ts_real &b) { return b < a; }
inline bool operator>(const ts_real &a, float b) {
  if (a[0] != b) return a[0] > b;
  if (a[1] != 0.0f) return a[1] > 0.0f;
  return a[2] > 0.0f;
}
inline bool operator>(float a, const ts_real &b) {
  if (a != b[0]) return a > b[0];
  if (b[1] != 0.0f) return 0.0f > b[1];
  return 0.0f > b[2];
}

inline bool operator<=(const ts_real &a, const ts_real &b) { return !(b < a); }
inline bool operator<=(const ts_real &a, float b)          { return !(a > b); }
inline bool operator<=(float a, const ts_real &b)          { return !(a > b); }
inline bool operator>=(const ts_real &a, const ts_real &b) { return !(a < b); }
inline bool operator>=(const ts_real &a, float b)          { return !(a < b); }
inline bool operator>=(float a, const ts_real &b)          { return !(a < b); }

/********** Misc **********/
inline ts_real abs(const ts_real &a)  { return (a[0] < 0.0f) ? -a : a; }
inline ts_real fabs(const ts_real &a) { return abs(a); }
inline ts_real sqr(const ts_real &a)  { return a * a; }

inline ts_real nint(const ts_real &a) {
  float c0 = qs::nint(a[0]);
  float c1 = 0.0f, c2 = 0.0f;
  if (c0 == a[0]) {
    c1 = qs::nint(a[1]);
    if (c1 == a[1]) c2 = qs::nint(a[2]);
    else if (std::fabs(c1 - a[1]) == 0.5f && a[2] < 0.0f) c1 -= 1.0f;
  } else if (std::fabs(c0 - a[0]) == 0.5f && a[1] < 0.0f) c0 -= 1.0f;
  qs::renorm3(c0, c1, c2);
  return ts_real(c0, c1, c2);
}

inline ts_real floor(const ts_real &a) {
  float c0 = std::floor(a[0]), c1 = 0.0f, c2 = 0.0f;
  if (c0 == a[0]) {
    c1 = std::floor(a[1]);
    if (c1 == a[1]) c2 = std::floor(a[2]);
  }
  qs::renorm3(c0, c1, c2);
  return ts_real(c0, c1, c2);
}
inline ts_real ceil(const ts_real &a) {
  float c0 = std::ceil(a[0]), c1 = 0.0f, c2 = 0.0f;
  if (c0 == a[0]) {
    c1 = std::ceil(a[1]);
    if (c1 == a[1]) c2 = std::ceil(a[2]);
  }
  qs::renorm3(c0, c1, c2);
  return ts_real(c0, c1, c2);
}
inline ts_real aint(const ts_real &a) { return (a[0] >= 0.0f) ? floor(a) : ceil(a); }

inline ds_real to_ds_real(const ts_real &a) { return ds_real(a[0], a[1]); }
inline float to_double(const ts_real &a) { return a[0]; }
inline float to_float (const ts_real &a) { return a[0]; }
inline int   to_int   (const ts_real &a) { return static_cast<int>(a[0]); }

inline ts_real ts_real::rand() { return tsrand(); }

inline ts_real max(const ts_real &a, const ts_real &b) { return (a > b) ? a : b; }
inline ts_real min(const ts_real &a, const ts_real &b) { return (a < b) ? a : b; }

#endif /* _QD_TS_REAL_H */
