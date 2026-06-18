/*
 * include/qd/fs_inline.h
 *
 * Basic error-free transformations (TwoSum, TwoProd, Veltkamp split, ...)
 * for IEEE-754 single precision (float).  These are the building blocks
 * used by the float-based ds_real / ts_real / qs_real classes.
 *
 * Created for the dtq package extension (2026).
 */
#ifndef _QD_FS_INLINE_H
#define _QD_FS_INLINE_H

#include <cmath>
#include <limits>
#include <qd/qd_config.h>

/* Splitter for Veltkamp split on a float (24-bit significand). */
#define _QS_SPLITTER     4097.0f               /* = 2^12 + 1 */
/* Threshold above which we have to scale before splitting to avoid overflow.
   float max ~ 3.4028e38 = 2^128, splitter ~ 2^12, so scale when |a| > 2^114. */
#define _QS_SPLIT_THRESH 2.0769187e+34f         /* = 2^114 */

namespace qs {

/* Quiet NaN / infinity for float. */
static const float _f_nan = std::numeric_limits<float>::quiet_NaN();
static const float _f_inf = std::numeric_limits<float>::infinity();

/* Computes fl(a+b) and err(a+b) assuming |a| >= |b|. */
inline float quick_two_sum(float a, float b, float &err) {
  float s = a + b;
  err = b - (s - a);
  return s;
}

/* Computes fl(a-b) and err(a-b) assuming |a| >= |b|. */
inline float quick_two_diff(float a, float b, float &err) {
  float s = a - b;
  err = (a - s) - b;
  return s;
}

/* Computes fl(a+b) and err(a+b). */
inline float two_sum(float a, float b, float &err) {
  float s = a + b;
  float bb = s - a;
  err = (a - (s - bb)) + (b - bb);
  return s;
}

/* Computes fl(a-b) and err(a-b). */
inline float two_diff(float a, float b, float &err) {
  float s = a - b;
  float bb = s - a;
  err = (a - (s - bb)) - (b + bb);
  return s;
}

/* Veltkamp split of a float into hi/lo halves of ~12 bits each. */
inline void split(float a, float &hi, float &lo) {
  float temp;
  if (a > _QS_SPLIT_THRESH || a < -_QS_SPLIT_THRESH) {
    a *= 1.220703125e-4f;             /* 2^-13 */
    temp = _QS_SPLITTER * a;
    hi = temp - (temp - a);
    lo = a - hi;
    hi *= 8192.0f;                    /* 2^13 */
    lo *= 8192.0f;
  } else {
    temp = _QS_SPLITTER * a;
    hi = temp - (temp - a);
    lo = a - hi;
  }
}

/* Computes fl(a*b) and err(a*b). */
inline float two_prod(float a, float b, float &err) {
#if defined(__FMA__) || defined(QD_FMS)
  /* If the platform offers a hardware FMA on float, use fmaf. */
  float p = a * b;
  err = std::fma(a, b, -p);
  return p;
#else
  float a_hi, a_lo, b_hi, b_lo;
  float p = a * b;
  split(a, a_hi, a_lo);
  split(b, b_hi, b_lo);
  err = ((a_hi * b_hi - p) + a_hi * b_lo + a_lo * b_hi) + a_lo * b_lo;
  return p;
#endif
}

/* Computes fl(a*a) and err(a*a). */
inline float two_sqr(float a, float &err) {
#if defined(__FMA__) || defined(QD_FMS)
  float p = a * a;
  err = std::fma(a, a, -p);
  return p;
#else
  float hi, lo;
  float q = a * a;
  split(a, hi, lo);
  err = ((hi * hi - q) + 2.0f * hi * lo) + lo * lo;
  return q;
#endif
}

/* Three-sum: (a + b + c) renormalized into (a, b, c). */
inline void three_sum(float &a, float &b, float &c) {
  float t1, t2, t3;
  t1 = qs::two_sum(a, b, t2);
  a  = qs::two_sum(c, t1, t3);
  b  = qs::two_sum(t2, t3, c);
}

/* Three-sum returning only the high two terms. */
inline void three_sum2(float &a, float &b, float &c) {
  float t1, t2, t3;
  t1 = qs::two_sum(a, b, t2);
  a  = qs::two_sum(c, t1, t3);
  b  = t2 + t3;
}

/* Nearest integer of a float. */
inline float nint(float d) {
  if (d == std::floor(d))
    return d;
  return std::floor(d + 0.5f);
}

/* Truncated integer of a float (toward zero). */
inline float aint(float d) {
  return (d >= 0.0f) ? std::floor(d) : std::ceil(d);
}

inline float sqr(float t) { return t * t; }

inline float to_float(float a) { return a; }

}

#endif /* _QD_FS_INLINE_H */
