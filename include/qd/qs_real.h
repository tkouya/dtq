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
inline qs_real operator/(const qs_real &a, const qs_real &b) {
  return qs_real::sloppy_div(a, b);
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
