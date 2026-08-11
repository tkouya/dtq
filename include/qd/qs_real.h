/*
 * include/qd/qs_real.h
 *
 * Quad-single precision (>= 96-bit significand) floating point arithmetic
 * package.  A qs_real is represented by an unevaluated sum of four IEEE
 * single-precision floats x[0] + x[1] + x[2] + x[3] in non-overlapping form,
 * giving roughly 28 decimal digits of precision (~4x float).
 *
 * Modeled directly on qd_real (LBNL).
 *
 * Created for the dtq package extension (2026).
 */
#ifndef _QD_QS_REAL_H
#define _QD_QS_REAL_H

#include <iostream>
#include <string>
#include <limits>
#include <qd/qd_config.h>
#include <qd/fs_inline.h>
#include <qd/ds_real.h>
#include <qd/ts_real.h>

/* Full ~96-bit precision holds only while the lowest limb (~|x|*2^-72)
   stays in the float normal range, i.e. |x| >= 2^-54; below that, precision
   degrades by about one bit per halving (representation limit of the
   format itself, independent of the algorithm). */
struct QD_API qs_real {
  float x[4];

  void renorm();
  void renorm(float &e);

  qs_real(float x0, float x1, float x2, float x3) {
    x[0] = x0; x[1] = x1; x[2] = x2; x[3] = x3;
  }
  explicit qs_real(const float *xx) {
    x[0] = xx[0]; x[1] = xx[1]; x[2] = xx[2]; x[3] = xx[3];
  }

  static const qs_real _2pi;
  static const qs_real _pi;
  static const qs_real _3pi4;
  static const qs_real _pi2;
  static const qs_real _pi4;
  static const qs_real _e;
  static const qs_real _log2;
  static const qs_real _log10;
  static const qs_real _nan;
  static const qs_real _inf;

  static const float _eps;
  static const float _min_normalized;
  static const qs_real _max;
  static const qs_real _safe_max;
  static const int _ndigits;

  qs_real() { x[0] = x[1] = x[2] = x[3] = 0.0f; }
  qs_real(const char *s);
  qs_real(const ds_real &a) { x[0] = a._hi(); x[1] = a._lo(); x[2] = x[3] = 0.0f; }
  qs_real(const ts_real &a) { x[0] = a[0]; x[1] = a[1]; x[2] = a[2]; x[3] = 0.0f; }
  qs_real(float d)  { x[0] = d;                                      x[1] = x[2] = x[3] = 0.0f; }
  qs_real(double d) { x[0] = static_cast<float>(d);                   x[1] = x[2] = x[3] = 0.0f; }
  qs_real(int i)    { x[0] = static_cast<float>(i);                   x[1] = x[2] = x[3] = 0.0f; }

  float operator[](int i) const { return x[i]; }
  float &operator[](int i) { return x[i]; }

  static void error(const char *msg);

  bool isnan() const {
    return QD_ISNAN(x[0]) || QD_ISNAN(x[1]) || QD_ISNAN(x[2]) || QD_ISNAN(x[3]);
  }
  bool isfinite() const { return QD_ISFINITE(x[0]); }
  bool isinf() const { return QD_ISINF(x[0]); }

  static qs_real ieee_add(const qs_real &a, const qs_real &b);
  static qs_real sloppy_add(const qs_real &a, const qs_real &b);

  qs_real &operator+=(float a);
  qs_real &operator+=(const ds_real &a);
  qs_real &operator+=(const qs_real &a);
  qs_real &operator-=(float a);
  qs_real &operator-=(const ds_real &a);
  qs_real &operator-=(const qs_real &a);

  static qs_real sloppy_mul(const qs_real &a, const qs_real &b);

  qs_real &operator*=(float a);
  qs_real &operator*=(const ds_real &a);
  qs_real &operator*=(const qs_real &a);

  static qs_real sloppy_div(const qs_real &a, const qs_real &b);
  static qs_real accurate_div(const qs_real &a, const qs_real &b);
  /* Division built on the branch-free quad-word FMA (qw_fma). */
  static qs_real fma_div(const qs_real &a, const qs_real &b);

  qs_real &operator/=(float a);
  qs_real &operator/=(const ds_real &a);
  qs_real &operator/=(const qs_real &a);

  qs_real operator^(int n) const;
  qs_real operator-() const { return qs_real(-x[0], -x[1], -x[2], -x[3]); }

  qs_real &operator=(float a)  { x[0] = a; x[1] = x[2] = x[3] = 0.0f; return *this; }
  qs_real &operator=(double a) { x[0] = static_cast<float>(a); x[1] = x[2] = x[3] = 0.0f; return *this; }
  qs_real &operator=(int a)    { x[0] = static_cast<float>(a); x[1] = x[2] = x[3] = 0.0f; return *this; }
  qs_real &operator=(const ds_real &a) { x[0] = a._hi(); x[1] = a._lo(); x[2] = x[3] = 0.0f; return *this; }
  qs_real &operator=(const ts_real &a) { x[0] = a[0]; x[1] = a[1]; x[2] = a[2]; x[3] = 0.0f; return *this; }
  qs_real &operator=(const char *s);

  bool is_zero()     const { return x[0] == 0.0f; }
  bool is_one()      const { return x[0] == 1.0f && x[1] == 0.0f && x[2] == 0.0f && x[3] == 0.0f; }
  bool is_positive() const { return x[0] >  0.0f; }
  bool is_negative() const { return x[0] <  0.0f; }

  static qs_real rand(void);

  std::string to_string(int precision = _ndigits) const;
};

namespace std {
  template <>
  class numeric_limits<qs_real> : public numeric_limits<float> {
  public:
    inline static float epsilon() { return qs_real::_eps; }
    inline static float min() { return qs_real::_min_normalized; }
    inline static qs_real max() { return qs_real::_max; }
    inline static qs_real safe_max() { return qs_real::_safe_max; }
    static const int digits = 96;
    static const int digits10 = 28;
  };
}

QD_API qs_real polyeval(const qs_real *c, int n, const qs_real &x);
QD_API qs_real polyroot(const qs_real *c, int n,
    const qs_real &x0, int max_iter = 64, float thresh = 0.0f);

QD_API qs_real qsrand(void);
QD_API qs_real sqrt(const qs_real &a);

QD_API inline bool isnan(const qs_real &a) { return a.isnan(); }
QD_API inline bool isfinite(const qs_real &a) { return a.isfinite(); }
QD_API inline bool isinf(const qs_real &a) { return a.isinf(); }

QD_API qs_real mul_pwr2(const qs_real &a, float b);

/*********** Branch-free quad-word FMA (T. Kouya) ************/
/* qw_fma(a, b, c) computes  a * b + c  as a single fused quad-word
   operation, without renormalizing a * b on its own. */
QD_API qs_real qw_fma(const qs_real &a, const qs_real &b, const qs_real &c);
QD_API qs_real qw_fma(const qs_real &a, float b, const qs_real &c);

/* Generic spelling; same operation. */
QD_API qs_real fma(const qs_real &a, const qs_real &b, const qs_real &c);
QD_API qs_real fma(const qs_real &a, float b, const qs_real &c);

/* Reference (pre-0.0.3) implementation, kept for benchmarking. */
QD_API qs_real sqrt_legacy(const qs_real &a);

QD_API qs_real operator+(const qs_real &a, const qs_real &b);
QD_API qs_real operator+(const qs_real &a, float b);
QD_API qs_real operator+(float a, const qs_real &b);

QD_API qs_real operator-(const qs_real &a, const qs_real &b);
QD_API qs_real operator-(const qs_real &a, float b);
QD_API qs_real operator-(float a, const qs_real &b);

QD_API qs_real operator*(const qs_real &a, const qs_real &b);
QD_API qs_real operator*(const qs_real &a, float b);
QD_API qs_real operator*(float a, const qs_real &b);

QD_API qs_real operator/(const qs_real &a, const qs_real &b);
QD_API qs_real operator/(const qs_real &a, float b);
QD_API qs_real operator/(float a, const qs_real &b);

QD_API qs_real sqr(const qs_real &a);
QD_API qs_real sqrt(const qs_real &a);
QD_API qs_real pow(const qs_real &a, int n);
QD_API qs_real pow(const qs_real &a, const qs_real &b);
QD_API qs_real npwr(const qs_real &a, int n);
QD_API qs_real nroot(const qs_real &a, int n);

ds_real to_ds_real(const qs_real &a);
ts_real to_ts_real(const qs_real &a);
QD_API float to_double(const qs_real &a);
QD_API float to_float (const qs_real &a);
QD_API int   to_int   (const qs_real &a);

QD_API bool operator==(const qs_real &a, const qs_real &b);
QD_API bool operator==(const qs_real &a, float b);
QD_API bool operator==(float a, const qs_real &b);

QD_API bool operator<(const qs_real &a, const qs_real &b);
QD_API bool operator<(const qs_real &a, float b);
QD_API bool operator<(float a, const qs_real &b);

QD_API bool operator>(const qs_real &a, const qs_real &b);
QD_API bool operator>(const qs_real &a, float b);
QD_API bool operator>(float a, const qs_real &b);

QD_API bool operator<=(const qs_real &a, const qs_real &b);
QD_API bool operator<=(const qs_real &a, float b);
QD_API bool operator<=(float a, const qs_real &b);

QD_API bool operator>=(const qs_real &a, const qs_real &b);
QD_API bool operator>=(const qs_real &a, float b);
QD_API bool operator>=(float a, const qs_real &b);

QD_API bool operator!=(const qs_real &a, const qs_real &b);
QD_API bool operator!=(const qs_real &a, float b);
QD_API bool operator!=(float a, const qs_real &b);

QD_API qs_real fabs(const qs_real &a);
QD_API qs_real abs(const qs_real &a);

QD_API qs_real ldexp(const qs_real &a, int n);
QD_API qs_real nint(const qs_real &a);
QD_API qs_real floor(const qs_real &a);
QD_API qs_real ceil(const qs_real &a);
QD_API qs_real aint(const qs_real &a);

QD_API qs_real exp(const qs_real &a);
QD_API qs_real log(const qs_real &a);
QD_API qs_real log10(const qs_real &a);

QD_API qs_real sin(const qs_real &a);
QD_API qs_real cos(const qs_real &a);
QD_API qs_real tan(const qs_real &a);
QD_API void    sincos(const qs_real &a, qs_real &s, qs_real &c);

QD_API std::ostream &operator<<(std::ostream &s, const qs_real &a);
QD_API std::istream &operator>>(std::istream &s, qs_real &a);

/********** Renormalization (in qs::) **********/
namespace qs {
inline void renorm4(float &c0, float &c1, float &c2, float &c3) {
  float s0, s1, s2 = 0.0f, s3 = 0.0f;
  if (QD_ISINF(c0)) return;
  s0 = qs::quick_two_sum(c2, c3, c3);
  s0 = qs::quick_two_sum(c1, s0, c2);
  c0 = qs::quick_two_sum(c0, s0, c1);
  s0 = c0; s1 = c1;
  if (s1 != 0.0f) {
    s1 = qs::quick_two_sum(s1, c2, s2);
    if (s2 != 0.0f) s2 = qs::quick_two_sum(s2, c3, s3);
    else            s1 = qs::quick_two_sum(s1, c3, s2);
  } else {
    s0 = qs::quick_two_sum(s0, c2, s1);
    if (s1 != 0.0f) s1 = qs::quick_two_sum(s1, c3, s2);
    else            s0 = qs::quick_two_sum(s0, c3, s1);
  }
  c0 = s0; c1 = s1; c2 = s2; c3 = s3;
}
inline void renorm5(float &c0, float &c1, float &c2, float &c3, float &c4) {
  float s0, s1, s2 = 0.0f, s3 = 0.0f;
  if (QD_ISINF(c0)) return;
  s0 = qs::quick_two_sum(c3, c4, c4);
  s0 = qs::quick_two_sum(c2, s0, c3);
  s0 = qs::quick_two_sum(c1, s0, c2);
  c0 = qs::quick_two_sum(c0, s0, c1);
  s0 = c0; s1 = c1;
  if (s1 != 0.0f) {
    s1 = qs::quick_two_sum(s1, c2, s2);
    if (s2 != 0.0f) {
      s2 = qs::quick_two_sum(s2, c3, s3);
      if (s3 != 0.0f) s3 += c4;
      else            s2 = qs::quick_two_sum(s2, c4, s3);
    } else {
      s1 = qs::quick_two_sum(s1, c3, s2);
      if (s2 != 0.0f) s2 = qs::quick_two_sum(s2, c4, s3);
      else            s1 = qs::quick_two_sum(s1, c4, s2);
    }
  } else {
    s0 = qs::quick_two_sum(s0, c2, s1);
    if (s1 != 0.0f) {
      s1 = qs::quick_two_sum(s1, c3, s2);
      if (s2 != 0.0f) s2 = qs::quick_two_sum(s2, c4, s3);
      else            s1 = qs::quick_two_sum(s1, c4, s2);
    } else {
      s0 = qs::quick_two_sum(s0, c3, s1);
      if (s1 != 0.0f) s1 = qs::quick_two_sum(s1, c4, s2);
      else            s0 = qs::quick_two_sum(s0, c4, s1);
    }
  }
  c0 = s0; c1 = s1; c2 = s2; c3 = s3;
}
}

inline void qs_real::renorm() {
  qs::renorm4(x[0], x[1], x[2], x[3]);
}
inline void qs_real::renorm(float &e) {
  qs::renorm5(x[0], x[1], x[2], x[3], e);
}

/********** Inline arithmetic **********/

inline qs_real operator+(const qs_real &a, float b) {
  float c0, c1, c2, c3, e;
  c0 = qs::two_sum(a[0], b, e);
  c1 = qs::two_sum(a[1], e, e);
  c2 = qs::two_sum(a[2], e, e);
  c3 = qs::two_sum(a[3], e, e);
  qs::renorm5(c0, c1, c2, c3, e);
  return qs_real(c0, c1, c2, c3);
}
inline qs_real operator+(float a, const qs_real &b) { return b + a; }

inline qs_real qs_real::sloppy_add(const qs_real &a, const qs_real &b) {
  float s0, s1, s2, s3;
  float t0, t1, t2, t3;
  s0 = qs::two_sum(a[0], b[0], t0);
  s1 = qs::two_sum(a[1], b[1], t1);
  s2 = qs::two_sum(a[2], b[2], t2);
  s3 = qs::two_sum(a[3], b[3], t3);
  s1 = qs::two_sum(s1, t0, t0);
  qs::three_sum(s2, t0, t1);
  qs::three_sum2(s3, t0, t2);
  t0 = t0 + t1 + t3;
  qs::renorm5(s0, s1, s2, s3, t0);
  return qs_real(s0, s1, s2, s3);
}

inline qs_real operator+(const qs_real &a, const qs_real &b) {
  return qs_real::sloppy_add(a, b);
}

inline qs_real operator-(const qs_real &a, const qs_real &b) { return a + (-b); }
inline qs_real operator-(const qs_real &a, float b)          { return a + (-b); }
inline qs_real operator-(float a, const qs_real &b)          { return a + (-b); }

inline qs_real operator*(const qs_real &a, float b) {
  float p0, p1, p2, p3;
  float q0, q1, q2;
  float s0, s1, s2, s3, s4;
  p0 = qs::two_prod(a[0], b, q0);
  p1 = qs::two_prod(a[1], b, q1);
  p2 = qs::two_prod(a[2], b, q2);
  p3 = a[3] * b;

  s0 = p0;
  s1 = qs::two_sum(q0, p1, s2);
  qs::three_sum(s2, q1, p2);
  qs::three_sum2(q1, q2, p3);
  s3 = q1;
  s4 = q2 + p2;
  qs::renorm5(s0, s1, s2, s3, s4);
  return qs_real(s0, s1, s2, s3);
}
inline qs_real operator*(float a, const qs_real &b) { return b * a; }

inline qs_real qs_real::sloppy_mul(const qs_real &a, const qs_real &b) {
  float p0, p1, p2, p3, p4, p5;
  float q0, q1, q2, q3, q4, q5;
  float t0, t1;
  float s0, s1, s2;

  p0 = qs::two_prod(a[0], b[0], q0);

  p1 = qs::two_prod(a[0], b[1], q1);
  p2 = qs::two_prod(a[1], b[0], q2);

  p3 = qs::two_prod(a[0], b[2], q3);
  p4 = qs::two_prod(a[1], b[1], q4);
  p5 = qs::two_prod(a[2], b[0], q5);

  /* Start accumulation */
  qs::three_sum(p1, p2, q0);

  /* Six-Three Sum  of  p2, q1, q2, p3, p4, p5 */
  qs::three_sum(p2, q1, q2);
  qs::three_sum(p3, p4, p5);
  /* compute (s0, s1, s2) = (p2, q1, q2) + (p3, p4, p5) */
  s0 = qs::two_sum(p2, p3, t0);
  s1 = qs::two_sum(q1, p4, t1);
  s2 = q2 + p5;
  s1 = qs::two_sum(s1, t0, t0);
  s2 += (t0 + t1);

  /* O(eps^3) order terms */
  s1 += a[0]*b[3] + a[1]*b[2] + a[2]*b[1] + a[3]*b[0] + q0 + q3 + q4 + q5;
  qs::renorm5(p0, p1, s0, s1, s2);
  return qs_real(p0, p1, s0, s1);
}
inline qs_real operator*(const qs_real &a, const qs_real &b) {
  return qs_real::sloppy_mul(a, b);
}

inline qs_real mul_pwr2(const qs_real &a, float b) {
  return qs_real(a[0]*b, a[1]*b, a[2]*b, a[3]*b);
}
inline qs_real ldexp(const qs_real &a, int n) {
  return qs_real(std::ldexp(a[0], n), std::ldexp(a[1], n),
                 std::ldexp(a[2], n), std::ldexp(a[3], n));
}

// 2026-08-04 T.Kouya
// Branch free algorithm: quad-word FMA,  z = a * b + c  (float limbs).
// Same network as the qd_real version; see include/qd/qd_inline.h.
/* QW-FMA  z = a * b + c   (176 flops)
   Machine-proved with FPANVerifier + z3 5.0.0 (ACS2026 formulation):
     error bound      |z-(ab+c)| <= 812 u^4 (|ab|+|c|)
     every FastTwoSum precondition  exp(x) >= exp(y)
     non-overlapping output         z0 |> z1 |> z2 |> z3   (strongly_dominates)
   The precision p is symbolic, so the proofs hold for binary32/64/128 alike.
   Normalization repeats a cascade over adjacent pairs; the pass count is the
   smallest for which the non-overlap is provable (DW 1 / TW 3 / QW 5). */
inline qs_real qw_fma(const qs_real &a, const qs_real &b, const qs_real &c) {
  float P00, E00, P01, E01, P10, E10, P02, E02, P11, E11, P20, E20;
  float P03, P12, P21, P30, D, B, r, Et;
  float A1, f1, f2, f3, f4;
  float A2, g1, g2, g3, g4, g5, g6, g7, g8, g9;
  float A3, t1, t2, t3, t4;
  float w0, w1, w2, w3;
  P00 = qs::two_prod(a[0], b[0], E00);
  P01 = qs::two_prod(a[0], b[1], E01);
  P10 = qs::two_prod(a[1], b[0], E10);
  P02 = qs::two_prod(a[0], b[2], E02);
  P11 = qs::two_prod(a[1], b[1], E11);
  P20 = qs::two_prod(a[2], b[0], E20);
  P03 = a[0] * b[3];
  P12 = a[1] * b[2];
  P21 = a[2] * b[1];
  P30 = a[3] * b[0];
  D   = (P03 + P30) + (P12 + P21);
  B   = qs::two_sum(P00, c[0], r);
  A1  = qs::two_sum(P01, P10, f1);
  A1  = qs::two_sum(A1, E00, f2);
  A1  = qs::two_sum(A1, c[1], f3);
  A1  = qs::two_sum(A1, r, f4);
  A2  = qs::two_sum(P02, P20, g1);
  A2  = qs::two_sum(A2, P11, g2);
  Et  = qs::two_sum(E01, E10, g4);
  A2  = qs::two_sum(A2, Et, g3);
  A2  = qs::two_sum(A2, c[2], g5);
  A2  = qs::two_sum(A2, f1, g6);
  A2  = qs::two_sum(A2, f2, g7);
  A2  = qs::two_sum(A2, f3, g8);
  A2  = qs::two_sum(A2, f4, g9);
  t1  = E02 + E20;
  t2  = E11 + D;
  t3  = (t1 + t2) + c[3];
  t1  = g1 + g2;
  t2  = g3 + g4;
  t1  = t1 + t2;
  t2  = g6 + g7;
  t4  = g8 + g9;
  t2  = t2 + t4;
  t1  = t1 + t2;
  t1  = t1 + g5;
  A3  = t3 + t1;
  float z0, z1, z2, z3;
  /* renormalization pass 1/5 */
  w0 = qs::quick_two_sum(B, A1, w1);
  w1 = qs::two_sum(w1, A2, w2);
  w2 = qs::two_sum(w2, A3, w3);
  /* renormalization pass 2/5 */
  w0 = qs::two_sum(w0, w1, w1);
  w1 = qs::two_sum(w1, w2, w2);
  w2 = qs::quick_two_sum(w2, w3, w3);
  /* renormalization pass 3/5 */
  w0 = qs::two_sum(w0, w1, w1);
  w1 = qs::quick_two_sum(w1, w2, w2);
  w2 = qs::quick_two_sum(w2, w3, w3);
  /* renormalization pass 4/5 */
  w0 = qs::quick_two_sum(w0, w1, w1);
  w1 = qs::quick_two_sum(w1, w2, w2);
  w2 = qs::quick_two_sum(w2, w3, w3);
  /* renormalization pass 5/5 */
  z0 = qs::quick_two_sum(w0, w1, w1);
  z1 = qs::quick_two_sum(w1, w2, w2);
  z2 = qs::quick_two_sum(w2, w3, z3);
  return qs_real(z0, z1, z2, z3);
}

/* QW-FMA  z = a * b + c   (176 flops, scalar multiplier)
   Machine-proved with FPANVerifier + z3 5.0.0 (ACS2026 formulation):
     error bound      |z-(ab+c)| <= 812 u^4 (|ab|+|c|)
     every FastTwoSum precondition  exp(x) >= exp(y)
     non-overlapping output         z0 |> z1 |> z2 |> z3   (strongly_dominates)
   The precision p is symbolic, so the proofs hold for binary32/64/128 alike.
   Normalization repeats a cascade over adjacent pairs; the pass count is the
   smallest for which the non-overlap is provable (DW 1 / TW 3 / QW 5). */

/* div/sqrt-safe variant (206 flops): the Newton iterations of division and
   square root receive residuals that are NOT non-overlapping expansions, so no
   FastTwoSum precondition can be claimed for any gate.  A FastTwoSum whose
   precondition fails does not even satisfy s+e=a+b -- it is no longer an EFT --
   so this variant uses TwoSum everywhere, with the same pass count as the
   standard variant.  A safe variant is never cheaper than the standard one. */
inline qs_real qw_fma_safe(const qs_real &a, float b, const qs_real &c) {
  float P00, E00, P01, E01, P10, E10, P02, E02, P11, E11, P20, E20;
  float P03, P12, P21, P30, D, B, r, Et;
  float A1, f1, f2, f3, f4;
  float A2, g1, g2, g3, g4, g5, g6, g7, g8, g9;
  float A3, t1, t2, t3, t4;
  float w0, w1, w2, w3;
  P00 = qs::two_prod(a[0], b, E00);
  P01 = 0.0; E01 = 0.0;
  P10 = qs::two_prod(a[1], b, E10);
  P02 = 0.0; E02 = 0.0;
  P11 = 0.0; E11 = 0.0;
  P20 = qs::two_prod(a[2], b, E20);
  P03 = 0.0;
  P12 = 0.0;
  P21 = 0.0;
  P30 = a[3] * b;
  D   = (P03 + P30) + (P12 + P21);
  B   = qs::two_sum(P00, c[0], r);
  A1  = qs::two_sum(P01, P10, f1);
  A1  = qs::two_sum(A1, E00, f2);
  A1  = qs::two_sum(A1, c[1], f3);
  A1  = qs::two_sum(A1, r, f4);
  A2  = qs::two_sum(P02, P20, g1);
  A2  = qs::two_sum(A2, P11, g2);
  Et  = qs::two_sum(E01, E10, g4);
  A2  = qs::two_sum(A2, Et, g3);
  A2  = qs::two_sum(A2, c[2], g5);
  A2  = qs::two_sum(A2, f1, g6);
  A2  = qs::two_sum(A2, f2, g7);
  A2  = qs::two_sum(A2, f3, g8);
  A2  = qs::two_sum(A2, f4, g9);
  t1  = E02 + E20;
  t2  = E11 + D;
  t3  = (t1 + t2) + c[3];
  t1  = g1 + g2;
  t2  = g3 + g4;
  t1  = t1 + t2;
  t2  = g6 + g7;
  t4  = g8 + g9;
  t2  = t2 + t4;
  t1  = t1 + t2;
  t1  = t1 + g5;
  A3  = t3 + t1;
  float z0, z1, z2, z3;
  /* renormalization pass 1/5 */
  w0 = qs::two_sum(B, A1, w1);
  w1 = qs::two_sum(w1, A2, w2);
  w2 = qs::two_sum(w2, A3, w3);
  /* renormalization pass 2/5 */
  w0 = qs::two_sum(w0, w1, w1);
  w1 = qs::two_sum(w1, w2, w2);
  w2 = qs::two_sum(w2, w3, w3);
  /* renormalization pass 3/5 */
  w0 = qs::two_sum(w0, w1, w1);
  w1 = qs::two_sum(w1, w2, w2);
  w2 = qs::two_sum(w2, w3, w3);
  /* renormalization pass 4/5 */
  w0 = qs::two_sum(w0, w1, w1);
  w1 = qs::two_sum(w1, w2, w2);
  w2 = qs::two_sum(w2, w3, w3);
  /* renormalization pass 5/5 */
  z0 = qs::two_sum(w0, w1, w1);
  z1 = qs::two_sum(w1, w2, w2);
  z2 = qs::two_sum(w2, w3, z3);
  return qs_real(z0, z1, z2, z3);
}

inline qs_real qw_fma(const qs_real &a, float b, const qs_real &c) {
  float P00, E00, P01, E01, P10, E10, P02, E02, P11, E11, P20, E20;
  float P03, P12, P21, P30, D, B, r, Et;
  float A1, f1, f2, f3, f4;
  float A2, g1, g2, g3, g4, g5, g6, g7, g8, g9;
  float A3, t1, t2, t3, t4;
  float w0, w1, w2, w3;
  P00 = qs::two_prod(a[0], b, E00);
  P01 = 0.0; E01 = 0.0;
  P10 = qs::two_prod(a[1], b, E10);
  P02 = 0.0; E02 = 0.0;
  P11 = 0.0; E11 = 0.0;
  P20 = qs::two_prod(a[2], b, E20);
  P03 = 0.0;
  P12 = 0.0;
  P21 = 0.0;
  P30 = a[3] * b;
  D   = (P03 + P30) + (P12 + P21);
  B   = qs::two_sum(P00, c[0], r);
  A1  = qs::two_sum(P01, P10, f1);
  A1  = qs::two_sum(A1, E00, f2);
  A1  = qs::two_sum(A1, c[1], f3);
  A1  = qs::two_sum(A1, r, f4);
  A2  = qs::two_sum(P02, P20, g1);
  A2  = qs::two_sum(A2, P11, g2);
  Et  = qs::two_sum(E01, E10, g4);
  A2  = qs::two_sum(A2, Et, g3);
  A2  = qs::two_sum(A2, c[2], g5);
  A2  = qs::two_sum(A2, f1, g6);
  A2  = qs::two_sum(A2, f2, g7);
  A2  = qs::two_sum(A2, f3, g8);
  A2  = qs::two_sum(A2, f4, g9);
  t1  = E02 + E20;
  t2  = E11 + D;
  t3  = (t1 + t2) + c[3];
  t1  = g1 + g2;
  t2  = g3 + g4;
  t1  = t1 + t2;
  t2  = g6 + g7;
  t4  = g8 + g9;
  t2  = t2 + t4;
  t1  = t1 + t2;
  t1  = t1 + g5;
  A3  = t3 + t1;
  float z0, z1, z2, z3;
  /* renormalization pass 1/5 */
  w0 = qs::quick_two_sum(B, A1, w1);
  w1 = qs::two_sum(w1, A2, w2);
  w2 = qs::two_sum(w2, A3, w3);
  /* renormalization pass 2/5 */
  w0 = qs::two_sum(w0, w1, w1);
  w1 = qs::two_sum(w1, w2, w2);
  w2 = qs::quick_two_sum(w2, w3, w3);
  /* renormalization pass 3/5 */
  w0 = qs::two_sum(w0, w1, w1);
  w1 = qs::quick_two_sum(w1, w2, w2);
  w2 = qs::quick_two_sum(w2, w3, w3);
  /* renormalization pass 4/5 */
  w0 = qs::quick_two_sum(w0, w1, w1);
  w1 = qs::quick_two_sum(w1, w2, w2);
  w2 = qs::quick_two_sum(w2, w3, w3);
  /* renormalization pass 5/5 */
  z0 = qs::quick_two_sum(w0, w1, w1);
  z1 = qs::quick_two_sum(w1, w2, w2);
  z2 = qs::quick_two_sum(w2, w3, z3);
  return qs_real(z0, z1, z2, z3);
}

inline qs_real fma(const qs_real &a, const qs_real &b, const qs_real &c) {
  return qw_fma(a, b, c);
}

inline qs_real fma(const qs_real &a, float b, const qs_real &c) {
  return qw_fma(a, b, c);
}

inline qs_real qs_real::sloppy_div(const qs_real &a, const qs_real &b) {
  float q0, q1, q2, q3;
  qs_real r;
  q0 = a[0] / b[0];
  r = a - b * q0;
  q1 = r[0] / b[0];
  r = r - b * q1;
  q2 = r[0] / b[0];
  r = r - b * q2;
  q3 = r[0] / b[0];
  qs::renorm4(q0, q1, q2, q3);
  return qs_real(q0, q1, q2, q3);
}
inline qs_real qs_real::accurate_div(const qs_real &a, const qs_real &b) {
  float q0, q1, q2, q3, q4;
  qs_real r;
  q0 = a[0] / b[0];
  r = a - b * q0;
  q1 = r[0] / b[0];
  r = r - b * q1;
  q2 = r[0] / b[0];
  r = r - b * q2;
  q3 = r[0] / b[0];
  r = r - b * q3;
  q4 = r[0] / b[0];
  qs::renorm5(q0, q1, q2, q3, q4);
  return qs_real(q0, q1, q2, q3);
}
/* Long division driven by the branch-free quad-word FMA.  Same correction
   sequence as sloppy_div (the qs_real default), but each residual
   r <- r - q * b  is one fused qw_fma. */
inline qs_real qs_real::fma_div(const qs_real &a, const qs_real &b) {
  float q0, q1, q2, q3;
  qs_real r;

  q0 = a[0] / b[0];
  r = qw_fma(b, -q0, a);
  q1 = r[0] / b[0];
  r = qw_fma(b, -q1, r);
  q2 = r[0] / b[0];
  r = qw_fma(b, -q2, r);
  q3 = r[0] / b[0];

  qs::renorm4(q0, q1, q2, q3);
  return qs_real(q0, q1, q2, q3);
}

inline qs_real operator/(const qs_real &a, const qs_real &b) {
#ifdef QD_NO_FMA_DIV
  return qs_real::sloppy_div(a, b);
#else
  return qs_real::fma_div(a, b);
#endif
}
inline qs_real operator/(const qs_real &a, float b) { return a / qs_real(b); }
inline qs_real operator/(float a, const qs_real &b) { return qs_real(a) / b; }

inline qs_real &qs_real::operator+=(float a)         { *this = *this + a; return *this; }
inline qs_real &qs_real::operator+=(const ds_real &a){ *this = *this + qs_real(a); return *this; }
inline qs_real &qs_real::operator+=(const qs_real &a){ *this = *this + a; return *this; }
inline qs_real &qs_real::operator-=(float a)         { *this = *this - a; return *this; }
inline qs_real &qs_real::operator-=(const ds_real &a){ *this = *this - qs_real(a); return *this; }
inline qs_real &qs_real::operator-=(const qs_real &a){ *this = *this - a; return *this; }
inline qs_real &qs_real::operator*=(float a)         { *this = *this * a; return *this; }
inline qs_real &qs_real::operator*=(const ds_real &a){ *this = *this * qs_real(a); return *this; }
inline qs_real &qs_real::operator*=(const qs_real &a){ *this = *this * a; return *this; }
inline qs_real &qs_real::operator/=(float a)         { *this = *this / a; return *this; }
inline qs_real &qs_real::operator/=(const ds_real &a){ *this = *this / qs_real(a); return *this; }
inline qs_real &qs_real::operator/=(const qs_real &a){ *this = *this / a; return *this; }

/********** Comparisons **********/
inline bool operator==(const qs_real &a, const qs_real &b) {
  return a[0] == b[0] && a[1] == b[1] && a[2] == b[2] && a[3] == b[3];
}
inline bool operator==(const qs_real &a, float b) {
  return a[0] == b && a[1] == 0.0f && a[2] == 0.0f && a[3] == 0.0f;
}
inline bool operator==(float a, const qs_real &b) { return b == a; }
inline bool operator!=(const qs_real &a, const qs_real &b) { return !(a == b); }
inline bool operator!=(const qs_real &a, float b)          { return !(a == b); }
inline bool operator!=(float a, const qs_real &b)          { return !(a == b); }

inline bool operator<(const qs_real &a, const qs_real &b) {
  if (a[0] != b[0]) return a[0] < b[0];
  if (a[1] != b[1]) return a[1] < b[1];
  if (a[2] != b[2]) return a[2] < b[2];
  return a[3] < b[3];
}
inline bool operator<(const qs_real &a, float b) {
  if (a[0] != b) return a[0] < b;
  if (a[1] != 0.0f) return a[1] < 0.0f;
  if (a[2] != 0.0f) return a[2] < 0.0f;
  return a[3] < 0.0f;
}
inline bool operator<(float a, const qs_real &b) {
  if (a != b[0]) return a < b[0];
  if (b[1] != 0.0f) return 0.0f < b[1];
  if (b[2] != 0.0f) return 0.0f < b[2];
  return 0.0f < b[3];
}

inline bool operator>(const qs_real &a, const qs_real &b) { return b < a; }
inline bool operator>(const qs_real &a, float b) {
  if (a[0] != b) return a[0] > b;
  if (a[1] != 0.0f) return a[1] > 0.0f;
  if (a[2] != 0.0f) return a[2] > 0.0f;
  return a[3] > 0.0f;
}
inline bool operator>(float a, const qs_real &b) {
  if (a != b[0]) return a > b[0];
  if (b[1] != 0.0f) return 0.0f > b[1];
  if (b[2] != 0.0f) return 0.0f > b[2];
  return 0.0f > b[3];
}

inline bool operator<=(const qs_real &a, const qs_real &b) { return !(b < a); }
inline bool operator<=(const qs_real &a, float b)          { return !(a > b); }
inline bool operator<=(float a, const qs_real &b)          { return !(a > b); }
inline bool operator>=(const qs_real &a, const qs_real &b) { return !(a < b); }
inline bool operator>=(const qs_real &a, float b)          { return !(a < b); }
inline bool operator>=(float a, const qs_real &b)          { return !(a < b); }

/********** Misc **********/
inline qs_real abs(const qs_real &a)  { return (a[0] < 0.0f) ? -a : a; }
inline qs_real fabs(const qs_real &a) { return abs(a); }
inline qs_real sqr(const qs_real &a)  { return a * a; }

inline qs_real nint(const qs_real &a) {
  float c0 = qs::nint(a[0]);
  float c1 = 0.0f, c2 = 0.0f, c3 = 0.0f;
  if (c0 == a[0]) {
    c1 = qs::nint(a[1]);
    if (c1 == a[1]) {
      c2 = qs::nint(a[2]);
      if (c2 == a[2]) c3 = qs::nint(a[3]);
    }
  }
  qs::renorm4(c0, c1, c2, c3);
  return qs_real(c0, c1, c2, c3);
}

inline qs_real floor(const qs_real &a) {
  float c0 = std::floor(a[0]), c1 = 0.0f, c2 = 0.0f, c3 = 0.0f;
  if (c0 == a[0]) {
    c1 = std::floor(a[1]);
    if (c1 == a[1]) {
      c2 = std::floor(a[2]);
      if (c2 == a[2]) c3 = std::floor(a[3]);
    }
  }
  qs::renorm4(c0, c1, c2, c3);
  return qs_real(c0, c1, c2, c3);
}
inline qs_real ceil(const qs_real &a) {
  float c0 = std::ceil(a[0]), c1 = 0.0f, c2 = 0.0f, c3 = 0.0f;
  if (c0 == a[0]) {
    c1 = std::ceil(a[1]);
    if (c1 == a[1]) {
      c2 = std::ceil(a[2]);
      if (c2 == a[2]) c3 = std::ceil(a[3]);
    }
  }
  qs::renorm4(c0, c1, c2, c3);
  return qs_real(c0, c1, c2, c3);
}
inline qs_real aint(const qs_real &a) { return (a[0] >= 0.0f) ? floor(a) : ceil(a); }

inline ds_real to_ds_real(const qs_real &a) { return ds_real(a[0], a[1]); }
inline ts_real to_ts_real(const qs_real &a) { return ts_real(a[0], a[1], a[2]); }
inline float to_double(const qs_real &a) { return a[0]; }
inline float to_float (const qs_real &a) { return a[0]; }
inline int   to_int   (const qs_real &a) { return static_cast<int>(a[0]); }

inline qs_real qs_real::rand() { return qsrand(); }

inline qs_real max(const qs_real &a, const qs_real &b) { return (a > b) ? a : b; }
inline qs_real min(const qs_real &a, const qs_real &b) { return (a < b) ? a : b; }

#endif /* _QD_QS_REAL_H */
