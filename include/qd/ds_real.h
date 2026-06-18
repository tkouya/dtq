/*
 * include/qd/ds_real.h
 *
 * Single-double precision (>= 48-bit significand) floating point
 * arithmetic package.  A ds_real is represented by an unevaluated sum
 * of two IEEE single-precision values x[0] + x[1] in non-overlapping
 * form, giving roughly 14 decimal digits of precision (~2x float).
 *
 * Modeled directly on dd_real (Y. Hida, X.S. Li, D.H. Bailey).
 *
 * Created for the dtq package extension (2026).
 */
#ifndef _QD_DS_REAL_H
#define _QD_DS_REAL_H

#include <cmath>
#include <iostream>
#include <string>
#include <limits>
#include <qd/qd_config.h>
#include <qd/fs_inline.h>

#ifdef isnan
#undef isnan
#endif
#ifdef isfinite
#undef isfinite
#endif
#ifdef isinf
#undef isinf
#endif
#ifdef max
#undef max
#endif
#ifdef min
#undef min
#endif

struct QD_API ds_real {
  float x[2];

  ds_real(float hi, float lo) { x[0] = hi; x[1] = lo; }
  ds_real() { x[0] = 0.0f; x[1] = 0.0f; }
  ds_real(float h) { x[0] = h; x[1] = 0.0f; }
  ds_real(double h) { x[0] = static_cast<float>(h); x[1] = 0.0f; }
  ds_real(int h) { x[0] = static_cast<float>(h); x[1] = 0.0f; }

  ds_real(const char *s);
  explicit ds_real(const float *d) { x[0] = d[0]; x[1] = d[1]; }

  static void error(const char *msg);

  float _hi() const { return x[0]; }
  float _lo() const { return x[1]; }

  static const ds_real _2pi;
  static const ds_real _pi;
  static const ds_real _3pi4;
  static const ds_real _pi2;
  static const ds_real _pi4;
  static const ds_real _e;
  static const ds_real _log2;
  static const ds_real _log10;
  static const ds_real _nan;
  static const ds_real _inf;

  static const float _eps;
  static const float _min_normalized;
  static const ds_real _max;
  static const ds_real _safe_max;
  static const int _ndigits;

  bool isnan() const { return QD_ISNAN(x[0]) || QD_ISNAN(x[1]); }
  bool isfinite() const { return QD_ISFINITE(x[0]); }
  bool isinf() const { return QD_ISINF(x[0]); }

  static ds_real add(float a, float b) {
    float s, e;
    s = qs::two_sum(a, b, e);
    return ds_real(s, e);
  }
  static ds_real sub(float a, float b) {
    float s, e;
    s = qs::two_diff(a, b, e);
    return ds_real(s, e);
  }
  static ds_real mul(float a, float b) {
    float p, e;
    p = qs::two_prod(a, b, e);
    return ds_real(p, e);
  }
  static ds_real div(float a, float b);
  static ds_real sqr(float a) {
    float p, e;
    p = qs::two_sqr(a, e);
    return ds_real(p, e);
  }
  static ds_real sqrt(float a);

  static ds_real ieee_add(const ds_real &a, const ds_real &b);
  static ds_real sloppy_add(const ds_real &a, const ds_real &b);
  static ds_real sloppy_div(const ds_real &a, const ds_real &b);
  static ds_real accurate_div(const ds_real &a, const ds_real &b);

  ds_real &operator+=(float a);
  ds_real &operator+=(const ds_real &a);
  ds_real &operator-=(float a);
  ds_real &operator-=(const ds_real &a);
  ds_real operator-() const { return ds_real(-x[0], -x[1]); }
  ds_real &operator*=(float a);
  ds_real &operator*=(const ds_real &a);
  ds_real &operator/=(float a);
  ds_real &operator/=(const ds_real &a);

  ds_real &operator=(float a) { x[0] = a; x[1] = 0.0f; return *this; }
  ds_real &operator=(double a) { x[0] = static_cast<float>(a); x[1] = 0.0f; return *this; }
  ds_real &operator=(int a)    { x[0] = static_cast<float>(a); x[1] = 0.0f; return *this; }
  ds_real &operator=(const char *s);

  ds_real operator^(int n) const;

  bool is_zero()     const { return x[0] == 0.0f; }
  bool is_one()      const { return x[0] == 1.0f && x[1] == 0.0f; }
  bool is_positive() const { return x[0] >  0.0f; }
  bool is_negative() const { return x[0] <  0.0f; }

  static ds_real rand(void);

  std::string to_string(int precision = _ndigits) const;
};

namespace std {
  template <>
  class numeric_limits<ds_real> : public numeric_limits<float> {
  public:
    inline static float epsilon() { return ds_real::_eps; }
    inline static ds_real max() { return ds_real::_max; }
    inline static ds_real safe_max() { return ds_real::_safe_max; }
    inline static float min() { return ds_real::_min_normalized; }
    static const int digits = 48;
    static const int digits10 = 14;
  };
}

QD_API ds_real dsrand(void);
QD_API ds_real sqrt(const ds_real &a);

QD_API ds_real polyeval(const ds_real *c, int n, const ds_real &x);
QD_API ds_real polyroot(const ds_real *c, int n,
    const ds_real &x0, int max_iter = 32, float thresh = 0.0f);

QD_API inline bool isnan(const ds_real &a) { return a.isnan(); }
QD_API inline bool isfinite(const ds_real &a) { return a.isfinite(); }
QD_API inline bool isinf(const ds_real &a) { return a.isinf(); }

QD_API ds_real mul_pwr2(const ds_real &a, float b);

QD_API ds_real operator+(const ds_real &a, float b);
QD_API ds_real operator+(float a, const ds_real &b);
QD_API ds_real operator+(const ds_real &a, const ds_real &b);

QD_API ds_real operator-(const ds_real &a, float b);
QD_API ds_real operator-(float a, const ds_real &b);
QD_API ds_real operator-(const ds_real &a, const ds_real &b);

QD_API ds_real operator*(const ds_real &a, float b);
QD_API ds_real operator*(float a, const ds_real &b);
QD_API ds_real operator*(const ds_real &a, const ds_real &b);

QD_API ds_real operator/(const ds_real &a, float b);
QD_API ds_real operator/(float a, const ds_real &b);
QD_API ds_real operator/(const ds_real &a, const ds_real &b);

QD_API ds_real inv(const ds_real &a);
QD_API ds_real rem(const ds_real &a, const ds_real &b);
QD_API ds_real drem(const ds_real &a, const ds_real &b);

QD_API ds_real pow(const ds_real &a, int n);
QD_API ds_real pow(const ds_real &a, const ds_real &b);
QD_API ds_real npwr(const ds_real &a, int n);
QD_API ds_real sqr(const ds_real &a);
QD_API ds_real sqrt(const ds_real &a);
QD_API ds_real nroot(const ds_real &a, int n);

QD_API bool operator==(const ds_real &a, float b);
QD_API bool operator==(float a, const ds_real &b);
QD_API bool operator==(const ds_real &a, const ds_real &b);

QD_API bool operator<=(const ds_real &a, float b);
QD_API bool operator<=(float a, const ds_real &b);
QD_API bool operator<=(const ds_real &a, const ds_real &b);

QD_API bool operator>=(const ds_real &a, float b);
QD_API bool operator>=(float a, const ds_real &b);
QD_API bool operator>=(const ds_real &a, const ds_real &b);

QD_API bool operator<(const ds_real &a, float b);
QD_API bool operator<(float a, const ds_real &b);
QD_API bool operator<(const ds_real &a, const ds_real &b);

QD_API bool operator>(const ds_real &a, float b);
QD_API bool operator>(float a, const ds_real &b);
QD_API bool operator>(const ds_real &a, const ds_real &b);

QD_API bool operator!=(const ds_real &a, float b);
QD_API bool operator!=(float a, const ds_real &b);
QD_API bool operator!=(const ds_real &a, const ds_real &b);

QD_API ds_real nint(const ds_real &a);
QD_API ds_real floor(const ds_real &a);
QD_API ds_real ceil(const ds_real &a);
QD_API ds_real aint(const ds_real &a);

QD_API ds_real fabs(const ds_real &a);
QD_API ds_real abs(const ds_real &a);

QD_API ds_real ldexp(const ds_real &a, int n);

QD_API ds_real exp(const ds_real &a);
QD_API ds_real log(const ds_real &a);
QD_API ds_real log10(const ds_real &a);

QD_API ds_real sin(const ds_real &a);
QD_API ds_real cos(const ds_real &a);
QD_API ds_real tan(const ds_real &a);
QD_API void    sincos(const ds_real &a, ds_real &s, ds_real &c);

QD_API float  to_double(const ds_real &a);   /* keep name to_double for template compatibility */
QD_API float  to_float(const ds_real &a);
QD_API int    to_int(const ds_real &a);

QD_API std::ostream &operator<<(std::ostream &s, const ds_real &a);
QD_API std::istream &operator>>(std::istream &s, ds_real &a);

/********** Inline arithmetic **********/
inline ds_real operator+(const ds_real &a, float b) {
  float s1, s2;
  s1 = qs::two_sum(a.x[0], b, s2);
  s2 += a.x[1];
  s1 = qs::quick_two_sum(s1, s2, s2);
  return ds_real(s1, s2);
}

inline ds_real operator+(float a, const ds_real &b) { return b + a; }

inline ds_real ds_real::sloppy_add(const ds_real &a, const ds_real &b) {
  float s, e;
  s = qs::two_sum(a.x[0], b.x[0], e);
  e += (a.x[1] + b.x[1]);
  s = qs::quick_two_sum(s, e, e);
  return ds_real(s, e);
}

inline ds_real ds_real::ieee_add(const ds_real &a, const ds_real &b) {
  float s1, s2, t1, t2;
  s1 = qs::two_sum(a.x[0], b.x[0], s2);
  t1 = qs::two_sum(a.x[1], b.x[1], t2);
  s2 += t1;
  s1 = qs::quick_two_sum(s1, s2, s2);
  s2 += t2;
  s1 = qs::quick_two_sum(s1, s2, s2);
  return ds_real(s1, s2);
}

inline ds_real operator+(const ds_real &a, const ds_real &b) {
  return ds_real::sloppy_add(a, b);
}

inline ds_real operator-(const ds_real &a, float b) {
  float s1, s2;
  s1 = qs::two_diff(a.x[0], b, s2);
  s2 += a.x[1];
  s1 = qs::quick_two_sum(s1, s2, s2);
  return ds_real(s1, s2);
}

inline ds_real operator-(float a, const ds_real &b) {
  float s1, s2;
  s1 = qs::two_diff(a, b.x[0], s2);
  s2 -= b.x[1];
  s1 = qs::quick_two_sum(s1, s2, s2);
  return ds_real(s1, s2);
}

inline ds_real operator-(const ds_real &a, const ds_real &b) {
  float s, e;
  s = qs::two_diff(a.x[0], b.x[0], e);
  e += a.x[1];
  e -= b.x[1];
  s = qs::quick_two_sum(s, e, e);
  return ds_real(s, e);
}

inline ds_real operator*(const ds_real &a, float b) {
  float p1, p2;
  p1 = qs::two_prod(a.x[0], b, p2);
  p2 += (a.x[1] * b);
  p1 = qs::quick_two_sum(p1, p2, p2);
  return ds_real(p1, p2);
}

inline ds_real operator*(float a, const ds_real &b) { return b * a; }

inline ds_real operator*(const ds_real &a, const ds_real &b) {
  float p1, p2;
  p1 = qs::two_prod(a.x[0], b.x[0], p2);
  p2 += (a.x[0] * b.x[1] + a.x[1] * b.x[0]);
  p1 = qs::quick_two_sum(p1, p2, p2);
  return ds_real(p1, p2);
}

inline ds_real mul_pwr2(const ds_real &a, float b) {
  return ds_real(a.x[0] * b, a.x[1] * b);
}

inline ds_real ldexp(const ds_real &a, int n) {
  return ds_real(std::ldexp(a.x[0], n), std::ldexp(a.x[1], n));
}

inline ds_real ds_real::sloppy_div(const ds_real &a, const ds_real &b) {
  float s1, s2;
  float q1, q2;
  ds_real r;
  q1 = a.x[0] / b.x[0];
  r = b * q1;
  s1 = qs::two_diff(a.x[0], r.x[0], s2);
  s2 -= r.x[1];
  s2 += a.x[1];
  q2 = (s1 + s2) / b.x[0];
  r.x[0] = qs::quick_two_sum(q1, q2, r.x[1]);
  return r;
}

inline ds_real ds_real::accurate_div(const ds_real &a, const ds_real &b) {
  float q1, q2, q3;
  ds_real r;
  q1 = a.x[0] / b.x[0];
  r = a - q1 * b;
  q2 = r.x[0] / b.x[0];
  r -= (q2 * b);
  q3 = r.x[0] / b.x[0];
  q1 = qs::quick_two_sum(q1, q2, q2);
  r = ds_real(q1, q2) + q3;
  return r;
}

inline ds_real operator/(const ds_real &a, float b) {
  float q1, q2;
  float p1, p2;
  float s, e;
  ds_real r;
  q1 = a.x[0] / b;
  p1 = qs::two_prod(q1, b, p2);
  s = qs::two_diff(a.x[0], p1, e);
  e += a.x[1];
  e -= p2;
  q2 = (s + e) / b;
  r.x[0] = qs::quick_two_sum(q1, q2, r.x[1]);
  return r;
}

inline ds_real operator/(const ds_real &a, const ds_real &b) {
  return ds_real::sloppy_div(a, b);
}

inline ds_real operator/(float a, const ds_real &b) {
  return ds_real(a) / b;
}

inline ds_real ds_real::div(float a, float b) {
  return ds_real(a) / ds_real(b);
}

inline ds_real inv(const ds_real &a) { return 1.0f / a; }

inline ds_real &ds_real::operator+=(float a) {
  float s1, s2;
  s1 = qs::two_sum(x[0], a, s2);
  s2 += x[1];
  x[0] = qs::quick_two_sum(s1, s2, x[1]);
  return *this;
}
inline ds_real &ds_real::operator+=(const ds_real &a) {
  float s, e;
  s = qs::two_sum(x[0], a.x[0], e);
  e += x[1]; e += a.x[1];
  x[0] = qs::quick_two_sum(s, e, x[1]);
  return *this;
}
inline ds_real &ds_real::operator-=(float a) {
  float s1, s2;
  s1 = qs::two_diff(x[0], a, s2);
  s2 += x[1];
  x[0] = qs::quick_two_sum(s1, s2, x[1]);
  return *this;
}
inline ds_real &ds_real::operator-=(const ds_real &a) {
  float s, e;
  s = qs::two_diff(x[0], a.x[0], e);
  e += x[1]; e -= a.x[1];
  x[0] = qs::quick_two_sum(s, e, x[1]);
  return *this;
}
inline ds_real &ds_real::operator*=(float a) {
  float p1, p2;
  p1 = qs::two_prod(x[0], a, p2);
  p2 += x[1] * a;
  x[0] = qs::quick_two_sum(p1, p2, x[1]);
  return *this;
}
inline ds_real &ds_real::operator*=(const ds_real &a) {
  float p1, p2;
  p1 = qs::two_prod(x[0], a.x[0], p2);
  p2 += a.x[1] * x[0];
  p2 += a.x[0] * x[1];
  x[0] = qs::quick_two_sum(p1, p2, x[1]);
  return *this;
}
inline ds_real &ds_real::operator/=(float a)         { *this = *this / a; return *this; }
inline ds_real &ds_real::operator/=(const ds_real &a) { *this = *this / a; return *this; }

/********** Comparisons **********/
inline bool operator==(const ds_real &a, float b)         { return a.x[0] == b && a.x[1] == 0.0f; }
inline bool operator==(float a, const ds_real &b)         { return a == b.x[0] && b.x[1] == 0.0f; }
inline bool operator==(const ds_real &a, const ds_real &b){ return a.x[0] == b.x[0] && a.x[1] == b.x[1]; }

inline bool operator<(const ds_real &a, float b)          { return a.x[0] < b || (a.x[0] == b && a.x[1] < 0.0f); }
inline bool operator<(float a, const ds_real &b)          { return a < b.x[0] || (a == b.x[0] && b.x[1] > 0.0f); }
inline bool operator<(const ds_real &a, const ds_real &b) { return a.x[0] < b.x[0] || (a.x[0] == b.x[0] && a.x[1] < b.x[1]); }

inline bool operator>(const ds_real &a, float b)          { return a.x[0] > b || (a.x[0] == b && a.x[1] > 0.0f); }
inline bool operator>(float a, const ds_real &b)          { return a > b.x[0] || (a == b.x[0] && b.x[1] < 0.0f); }
inline bool operator>(const ds_real &a, const ds_real &b) { return a.x[0] > b.x[0] || (a.x[0] == b.x[0] && a.x[1] > b.x[1]); }

inline bool operator<=(const ds_real &a, float b)          { return a.x[0] < b || (a.x[0] == b && a.x[1] <= 0.0f); }
inline bool operator<=(float a, const ds_real &b)          { return b >= a; }
inline bool operator<=(const ds_real &a, const ds_real &b) { return a.x[0] < b.x[0] || (a.x[0] == b.x[0] && a.x[1] <= b.x[1]); }

inline bool operator>=(const ds_real &a, float b)          { return a.x[0] > b || (a.x[0] == b && a.x[1] >= 0.0f); }
inline bool operator>=(float a, const ds_real &b)          { return b <= a; }
inline bool operator>=(const ds_real &a, const ds_real &b) { return a.x[0] > b.x[0] || (a.x[0] == b.x[0] && a.x[1] >= b.x[1]); }

inline bool operator!=(const ds_real &a, float b)          { return !(a == b); }
inline bool operator!=(float a, const ds_real &b)          { return !(a == b); }
inline bool operator!=(const ds_real &a, const ds_real &b) { return !(a == b); }

/********** Misc **********/
inline ds_real abs(const ds_real &a) { return (a.x[0] < 0.0f) ? -a : a; }
inline ds_real fabs(const ds_real &a) { return abs(a); }
inline ds_real sqr(const ds_real &a) {
  float p1, p2, s1, s2;
  p1 = qs::two_sqr(a.x[0], p2);
  p2 += 2.0f * a.x[0] * a.x[1];
  p2 += a.x[1] * a.x[1];
  s1 = qs::quick_two_sum(p1, p2, s2);
  return ds_real(s1, s2);
}

inline ds_real nint(const ds_real &a) {
  float hi = qs::nint(a.x[0]);
  float lo;
  if (hi == a.x[0]) {
    lo = qs::nint(a.x[1]);
    hi = qs::quick_two_sum(hi, lo, lo);
  } else {
    lo = 0.0f;
    if (std::fabs(hi - a.x[0]) == 0.5f && a.x[1] < 0.0f) hi -= 1.0f;
  }
  return ds_real(hi, lo);
}

inline ds_real floor(const ds_real &a) {
  float hi = std::floor(a.x[0]);
  float lo = 0.0f;
  if (hi == a.x[0]) {
    lo = std::floor(a.x[1]);
    hi = qs::quick_two_sum(hi, lo, lo);
  }
  return ds_real(hi, lo);
}

inline ds_real ceil(const ds_real &a) {
  float hi = std::ceil(a.x[0]);
  float lo = 0.0f;
  if (hi == a.x[0]) {
    lo = std::ceil(a.x[1]);
    hi = qs::quick_two_sum(hi, lo, lo);
  }
  return ds_real(hi, lo);
}

inline ds_real aint(const ds_real &a) {
  return (a.x[0] >= 0.0f) ? floor(a) : ceil(a);
}

inline float to_double(const ds_real &a) { return a.x[0]; }
inline float to_float (const ds_real &a) { return a.x[0]; }
inline int   to_int   (const ds_real &a) { return static_cast<int>(a.x[0]); }

inline ds_real ds_real::rand() { return dsrand(); }

inline ds_real drem(const ds_real &a, const ds_real &b) {
  ds_real n = nint(a / b);
  return a - n * b;
}

#endif /* _QD_DS_REAL_H */
