/*
 * include/dd_inline.h
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2000-2001
 *
 * Contains small functions (suitable for inlining) in the double-double
 * arithmetic package.
 */
#ifndef _QD_DD_INLINE_H
#define _QD_DD_INLINE_H

#include <cmath>
#include <qd/inline.h>

#ifndef QD_INLINE
#define inline
#endif


/*********** Additions ************/
/* double-double = double + double */
inline dd_real dd_real::add(double a, double b) {
  double s, e;
  s = qd::two_sum(a, b, e);
  return dd_real(s, e);
}

/* double-double + double */
inline dd_real operator+(const dd_real &a, double b) {
  double s1, s2;
  s1 = qd::two_sum(a.x[0], b, s2);
  s2 += a.x[1];
  s1 = qd::quick_two_sum(s1, s2, s2);
  return dd_real(s1, s2);
}

/* double-double + double-double */
inline dd_real dd_real::ieee_add(const dd_real &a, const dd_real &b) {
  /* This one satisfies IEEE style error bound, 
     due to K. Briggs and W. Kahan.                   */
  double s1, s2, t1, t2;

  s1 = qd::two_sum(a.x[0], b.x[0], s2);
  t1 = qd::two_sum(a.x[1], b.x[1], t2);
  s2 += t1;
  s1 = qd::quick_two_sum(s1, s2, s2);
  s2 += t2;
  s1 = qd::quick_two_sum(s1, s2, s2);
  return dd_real(s1, s2);
}

inline dd_real dd_real::sloppy_add(const dd_real &a, const dd_real &b) {
  /* This is the less accurate version ... obeys Cray-style
     error bound. */
  double s, e;

  s = qd::two_sum(a.x[0], b.x[0], e);
  e += (a.x[1] + b.x[1]);
  s = qd::quick_two_sum(s, e, e);
  return dd_real(s, e);
}

inline dd_real operator+(const dd_real &a, const dd_real &b) {
#ifndef QD_IEEE_ADD
  return dd_real::sloppy_add(a, b);
#else
  return dd_real::ieee_add(a, b);
#endif
}

/* double + double-double */
inline dd_real operator+(double a, const dd_real &b) {
  return (b + a);
}


/*********** Self-Additions ************/
/* double-double += double */
inline dd_real &dd_real::operator+=(double a) {
  double s1, s2;
  s1 = qd::two_sum(x[0], a, s2);
  s2 += x[1];
  x[0] = qd::quick_two_sum(s1, s2, x[1]);
  return *this;
}

/* double-double += double-double */
inline dd_real &dd_real::operator+=(const dd_real &a) {
#ifndef QD_IEEE_ADD
  double s, e;
  s = qd::two_sum(x[0], a.x[0], e);
  e += x[1];
  e += a.x[1];
  x[0] = qd::quick_two_sum(s, e, x[1]);
  return *this;
#else
  double s1, s2, t1, t2;
  s1 = qd::two_sum(x[0], a.x[0], s2);
  t1 = qd::two_sum(x[1], a.x[1], t2);
  s2 += t1;
  s1 = qd::quick_two_sum(s1, s2, s2);
  s2 += t2;
  x[0] = qd::quick_two_sum(s1, s2, x[1]);
  return *this;
#endif
}

/*********** Subtractions ************/
/* double-double = double - double */
inline dd_real dd_real::sub(double a, double b) {
  double s, e;
  s = qd::two_diff(a, b, e);
  return dd_real(s, e);
}

/* double-double - double */
inline dd_real operator-(const dd_real &a, double b) {
  double s1, s2;
  s1 = qd::two_diff(a.x[0], b, s2);
  s2 += a.x[1];
  s1 = qd::quick_two_sum(s1, s2, s2);
  return dd_real(s1, s2);
}

/* double-double - double-double */
inline dd_real operator-(const dd_real &a, const dd_real &b) {
#ifndef QD_IEEE_ADD
  double s, e;
  s = qd::two_diff(a.x[0], b.x[0], e);
  e += a.x[1];
  e -= b.x[1];
  s = qd::quick_two_sum(s, e, e);
  return dd_real(s, e);
#else
  double s1, s2, t1, t2;
  s1 = qd::two_diff(a.x[0], b.x[0], s2);
  t1 = qd::two_diff(a.x[1], b.x[1], t2);
  s2 += t1;
  s1 = qd::quick_two_sum(s1, s2, s2);
  s2 += t2;
  s1 = qd::quick_two_sum(s1, s2, s2);
  return dd_real(s1, s2);
#endif
}

/* double - double-double */
inline dd_real operator-(double a, const dd_real &b) {
  double s1, s2;
  s1 = qd::two_diff(a, b.x[0], s2);
  s2 -= b.x[1];
  s1 = qd::quick_two_sum(s1, s2, s2);
  return dd_real(s1, s2);
}

/*********** Self-Subtractions ************/
/* double-double -= double */
inline dd_real &dd_real::operator-=(double a) {
  double s1, s2;
  s1 = qd::two_diff(x[0], a, s2);
  s2 += x[1];
  x[0] = qd::quick_two_sum(s1, s2, x[1]);
  return *this;
}

/* double-double -= double-double */
inline dd_real &dd_real::operator-=(const dd_real &a) {
#ifndef QD_IEEE_ADD
  double s, e;
  s = qd::two_diff(x[0], a.x[0], e);
  e += x[1];
  e -= a.x[1];
  x[0] = qd::quick_two_sum(s, e, x[1]);
  return *this;
#else
  double s1, s2, t1, t2;
  s1 = qd::two_diff(x[0], a.x[0], s2);
  t1 = qd::two_diff(x[1], a.x[1], t2);
  s2 += t1;
  s1 = qd::quick_two_sum(s1, s2, s2);
  s2 += t2;
  x[0] = qd::quick_two_sum(s1, s2, x[1]);
  return *this;
#endif
}

/*********** Unary Minus ***********/
inline dd_real dd_real::operator-() const {
  return dd_real(-x[0], -x[1]);
}

/*********** Multiplications ************/
/* double-double = double * double */
inline dd_real dd_real::mul(double a, double b) {
  double p, e;
  p = qd::two_prod(a, b, e);
  return dd_real(p, e);
}

/* double-double * (2.0 ^ exp) */
inline dd_real ldexp(const dd_real &a, int exp) {
  return dd_real(std::ldexp(a.x[0], exp), std::ldexp(a.x[1], exp));
}

/* double-double * double,  where double is a power of 2. */
inline dd_real mul_pwr2(const dd_real &a, double b) {
  return dd_real(a.x[0] * b, a.x[1] * b);
}

/* double-double * double */
inline dd_real operator*(const dd_real &a, double b) {
  double p1, p2;

  p1 = qd::two_prod(a.x[0], b, p2);
  p2 += (a.x[1] * b);
  p1 = qd::quick_two_sum(p1, p2, p2);
  return dd_real(p1, p2);
}

/* double-double * double-double */
inline dd_real operator*(const dd_real &a, const dd_real &b) {
  double p1, p2;

  p1 = qd::two_prod(a.x[0], b.x[0], p2);
  p2 += (a.x[0] * b.x[1] + a.x[1] * b.x[0]);
  p1 = qd::quick_two_sum(p1, p2, p2);
  return dd_real(p1, p2);
}

/* double * double-double */
inline dd_real operator*(double a, const dd_real &b) {
  return (b * a);
}

/*********** Self-Multiplications ************/
/* double-double *= double */
inline dd_real &dd_real::operator*=(double a) {
  double p1, p2;
  p1 = qd::two_prod(x[0], a, p2);
  p2 += x[1] * a;
  x[0] = qd::quick_two_sum(p1, p2, x[1]);
  return *this;
}

/* double-double *= double-double */
inline dd_real &dd_real::operator*=(const dd_real &a) {
  double p1, p2;
  p1 = qd::two_prod(x[0], a.x[0], p2);
  p2 += a.x[1] * x[0];
  p2 += a.x[0] * x[1];
  x[0] = qd::quick_two_sum(p1, p2, x[1]);
  return *this;
}

/*********** Branch-free fused multiply-add ************/
// 2026-08-04 T.Kouya
// Branch free algorithm: double-word FMA,  z = a * b + c.
//
// The three inputs are split over two "levels" of magnitude (level 0 is
// the leading word, level 1 the trailing word).  Every term of the exact
// product a*b and of c is dropped into the level it belongs to, the two
// levels are accumulated with two_sum, and the result is renormalized with
// a single quick_two_sum.  No branch, and no intermediate renormalization
// of a*b.
/* DW-FMA  z = a * b + c   (17 flops)
   Machine-proved with FPANVerifier + z3 5.0.0 (ACS2026 formulation):
     error bound      |z-(ab+c)| <= 34 u^2 (|ab|+|c|)
     every FastTwoSum precondition  exp(x) >= exp(y)
     non-overlapping output         z0 |> z1   (strongly_dominates)
   The precision p is symbolic, so the proofs hold for binary32/64/128 alike. */
inline dd_real dw_fma(const dd_real &a, const dd_real &b, const dd_real &c) {
  double P00, E00, P01, P10, l, v, w, s, t, tp;
  P00 = qd::two_prod(a.x[0], b.x[0], E00);
  P01 = a.x[0] * b.x[1];
  P10 = a.x[1] * b.x[0];
  l   = P01 + P10;
  v   = E00 + c.x[1];
  w   = v + l;
  s   = qd::two_sum(P00, c.x[0], t);
  tp  = t + w;
  double z0, z1;
  z0 = qd::quick_two_sum(s, tp, z1);
  return dd_real(z0, z1);
}

/* DW-FMA  z = a * b + c   (17 flops, scalar multiplier)
   Machine-proved with FPANVerifier + z3 5.0.0 (ACS2026 formulation):
     error bound      |z-(ab+c)| <= 34 u^2 (|ab|+|c|)
     every FastTwoSum precondition  exp(x) >= exp(y)
     non-overlapping output         z0 |> z1   (strongly_dominates)
   The precision p is symbolic, so the proofs hold for binary32/64/128 alike. */

/* div/sqrt-safe variant: the Newton iterations of division and square root
   receive residuals that are NOT non-overlapping expansions, so no FastTwoSum
   precondition can be claimed for any gate.  A FastTwoSum whose precondition
   fails does not even satisfy s+e=a+b, so this variant uses TwoSum everywhere
   (20 flops); a safe variant is never cheaper than the standard one. */
inline dd_real dw_fma_safe(const dd_real &a, double b, const dd_real &c) {
  double P00, E00, P01, P10, l, v, w, s, t, tp;
  P00 = qd::two_prod(a.x[0], b, E00);
  P01 = 0.0;
  P10 = a.x[1] * b;
  l   = P01 + P10;
  v   = E00 + c.x[1];
  w   = v + l;
  s   = qd::two_sum(P00, c.x[0], t);
  tp  = t + w;
  double z0, z1;
  z0 = qd::two_sum(s, tp, z1);
  return dd_real(z0, z1);
}

inline dd_real dw_fma(const dd_real &a, double b, const dd_real &c) {
  double P00, E00, P01, P10, l, v, w, s, t, tp;
  P00 = qd::two_prod(a.x[0], b, E00);
  P01 = 0.0;
  P10 = a.x[1] * b;
  l   = P01 + P10;
  v   = E00 + c.x[1];
  w   = v + l;
  s   = qd::two_sum(P00, c.x[0], t);
  tp  = t + w;
  double z0, z1;
  z0 = qd::quick_two_sum(s, tp, z1);
  return dd_real(z0, z1);
}

/* double-word FMA with two plain double multiplicands:  a * b + c. */
inline dd_real dw_fma(double a, double b, const dd_real &c) {
  double a0, b0, d0;
  double a1, b1;
  double a2, b2;

  a0 = qd::two_prod(a, b, b0);
  a1 = qd::two_sum(a0, c.x[0], b1);
  d0 = (b0 + c.x[1]) + b1;
  a2 = qd::quick_two_sum(a1, d0, b2);

  return dd_real(a2, b2);
}

inline dd_real fma(const dd_real &a, const dd_real &b, const dd_real &c) {
  return dw_fma(a, b, c);
}

inline dd_real fma(const dd_real &a, double b, const dd_real &c) {
  return dw_fma(a, b, c);
}

/*********** Divisions ************/
inline dd_real dd_real::div(double a, double b) {
  double q1, q2;
  double p1, p2;
  double s, e;

  q1 = a / b;

  /* Compute  a - q1 * b */
  p1 = qd::two_prod(q1, b, p2);
  s = qd::two_diff(a, p1, e);
  e -= p2;

  /* get next approximation */
  q2 = (s + e) / b;

  s = qd::quick_two_sum(q1, q2, e);

  return dd_real(s, e);
}

/* double-double / double */
inline dd_real operator/(const dd_real &a, double b) {

  double q1, q2;
  double p1, p2;
  double s, e;
  dd_real r;
  
  q1 = a.x[0] / b;   /* approximate quotient. */

  /* Compute  this - q1 * d */
  p1 = qd::two_prod(q1, b, p2);
  s = qd::two_diff(a.x[0], p1, e);
  e += a.x[1];
  e -= p2;
  
  /* get next approximation. */
  q2 = (s + e) / b;

  /* renormalize */
  r.x[0] = qd::quick_two_sum(q1, q2, r.x[1]);

  return r;
}

inline dd_real dd_real::sloppy_div(const dd_real &a, const dd_real &b) {
  double s1, s2;
  double q1, q2;
  dd_real r;

  q1 = a.x[0] / b.x[0];  /* approximate quotient */

  /* compute  this - q1 * dd */
  r = b * q1;
  s1 = qd::two_diff(a.x[0], r.x[0], s2);
  s2 -= r.x[1];
  s2 += a.x[1];

  /* get next approximation */
  q2 = (s1 + s2) / b.x[0];

  /* renormalize */
  r.x[0] = qd::quick_two_sum(q1, q2, r.x[1]);
  return r;
}

inline dd_real dd_real::accurate_div(const dd_real &a, const dd_real &b) {
  double q1, q2, q3;
  dd_real r;

  q1 = a.x[0] / b.x[0];  /* approximate quotient */

  r = a - q1 * b;
  
  q2 = r.x[0] / b.x[0];
  r -= (q2 * b);

  q3 = r.x[0] / b.x[0];

  q1 = qd::quick_two_sum(q1, q2, q2);
  r = dd_real(q1, q2) + q3;
  return r;
}

/* Long division driven by the branch-free double-word FMA.  Same
   correction sequence as accurate_div, but each residual
   r <- r - q * b  is one fused dw_fma instead of a multiply followed by
   a subtraction, which removes one renormalization per step. */
inline dd_real dd_real::fma_div(const dd_real &a, const dd_real &b) {
  double q1, q2, q3;
  dd_real r;

  q1 = a.x[0] / b.x[0];      /* approximate quotient */

  r = dw_fma(b, -q1, a);     /* r = a - q1 * b */

  q2 = r.x[0] / b.x[0];
  r = dw_fma(b, -q2, r);     /* r = r - q2 * b */

  q3 = r.x[0] / b.x[0];

  q1 = qd::quick_two_sum(q1, q2, q2);
  return dd_real(q1, q2) + q3;
}

/* double-double / double-double */
inline dd_real operator/(const dd_real &a, const dd_real &b) {
#ifdef QD_SLOPPY_DIV
  return dd_real::sloppy_div(a, b);
#elif defined(QD_NO_FMA_DIV)
  return dd_real::accurate_div(a, b);
#else
  return dd_real::fma_div(a, b);
#endif
}

/* double / double-double */
inline dd_real operator/(double a, const dd_real &b) {
  return dd_real(a) / b;
}

inline dd_real inv(const dd_real &a) {
  return 1.0 / a;
}

/*********** Self-Divisions ************/
/* double-double /= double */
inline dd_real &dd_real::operator/=(double a) {
  *this = *this / a;
  return *this;
}

/* double-double /= double-double */
inline dd_real &dd_real::operator/=(const dd_real &a) {
  *this = *this / a;
  return *this;
}

/********** Remainder **********/
inline dd_real drem(const dd_real &a, const dd_real &b) {
  dd_real n = nint(a / b);
  return (a - n * b);
}

inline dd_real divrem(const dd_real &a, const dd_real &b, dd_real &r) {
  dd_real n = nint(a / b);
  r = a - n * b;
  return n;
}

/*********** Squaring **********/
inline dd_real sqr(const dd_real &a) {
  double p1, p2;
  double s1, s2;
  p1 = qd::two_sqr(a.x[0], p2);
  p2 += 2.0 * a.x[0] * a.x[1];
  p2 += a.x[1] * a.x[1];
  s1 = qd::quick_two_sum(p1, p2, s2);
  return dd_real(s1, s2);
}

inline dd_real dd_real::sqr(double a) {
  double p1, p2;
  p1 = qd::two_sqr(a, p2);
  return dd_real(p1, p2);
}


/********** Exponentiation **********/
inline dd_real dd_real::operator^(int n) {
  return npwr(*this, n);
}


/*********** Assignments ************/
/* double-double = double */
inline dd_real &dd_real::operator=(double a) {
  x[0] = a;
  x[1] = 0.0;
  return *this;
}

/*********** Equality Comparisons ************/
/* double-double == double */
inline bool operator==(const dd_real &a, double b) {
  return (a.x[0] == b && a.x[1] == 0.0);
}

/* double-double == double-double */
inline bool operator==(const dd_real &a, const dd_real &b) {
  return (a.x[0] == b.x[0] && a.x[1] == b.x[1]);
}

/* double == double-double */
inline bool operator==(double a, const dd_real &b) {
  return (a == b.x[0] && b.x[1] == 0.0);
}

/*********** Greater-Than Comparisons ************/
/* double-double > double */
inline bool operator>(const dd_real &a, double b) {
  return (a.x[0] > b || (a.x[0] == b && a.x[1] > 0.0));
}

/* double-double > double-double */
inline bool operator>(const dd_real &a, const dd_real &b) {
  return (a.x[0] > b.x[0] || (a.x[0] == b.x[0] && a.x[1] > b.x[1]));
}

/* double > double-double */
inline bool operator>(double a, const dd_real &b) {
  return (a > b.x[0] || (a == b.x[0] && b.x[1] < 0.0));
}

/*********** Less-Than Comparisons ************/
/* double-double < double */
inline bool operator<(const dd_real &a, double b) {
  return (a.x[0] < b || (a.x[0] == b && a.x[1] < 0.0));
}

/* double-double < double-double */
inline bool operator<(const dd_real &a, const dd_real &b) {
  return (a.x[0] < b.x[0] || (a.x[0] == b.x[0] && a.x[1] < b.x[1]));
}

/* double < double-double */
inline bool operator<(double a, const dd_real &b) {
  return (a < b.x[0] || (a == b.x[0] && b.x[1] > 0.0));
}

/*********** Greater-Than-Or-Equal-To Comparisons ************/
/* double-double >= double */
inline bool operator>=(const dd_real &a, double b) {
  return (a.x[0] > b || (a.x[0] == b && a.x[1] >= 0.0));
}

/* double-double >= double-double */
inline bool operator>=(const dd_real &a, const dd_real &b) {
  return (a.x[0] > b.x[0] || (a.x[0] == b.x[0] && a.x[1] >= b.x[1]));
}

/* double >= double-double */
inline bool operator>=(double a, const dd_real &b) {
  return (b <= a);
}

/*********** Less-Than-Or-Equal-To Comparisons ************/
/* double-double <= double */
inline bool operator<=(const dd_real &a, double b) {
  return (a.x[0] < b || (a.x[0] == b && a.x[1] <= 0.0));
}

/* double-double <= double-double */
inline bool operator<=(const dd_real &a, const dd_real &b) {
  return (a.x[0] < b.x[0] || (a.x[0] == b.x[0] && a.x[1] <= b.x[1]));
}

/* double <= double-double */
inline bool operator<=(double a, const dd_real &b) {
  return (b >= a);
}

/*********** Not-Equal-To Comparisons ************/
/* double-double != double */
inline bool operator!=(const dd_real &a, double b) {
  return (a.x[0] != b || a.x[1] != 0.0);
}

/* double-double != double-double */
inline bool operator!=(const dd_real &a, const dd_real &b) {
  return (a.x[0] != b.x[0] || a.x[1] != b.x[1]);
}

/* double != double-double */
inline bool operator!=(double a, const dd_real &b) {
  return (a != b.x[0] || b.x[1] != 0.0);
}

/*********** Micellaneous ************/
/*  this == 0 */
inline bool dd_real::is_zero() const {
  return (x[0] == 0.0);
}

/*  this == 1 */
inline bool dd_real::is_one() const {
  return (x[0] == 1.0 && x[1] == 0.0);
}

/*  this > 0 */
inline bool dd_real::is_positive() const {
  return (x[0] > 0.0);
}

/* this < 0 */
inline bool dd_real::is_negative() const {
  return (x[0] < 0.0);
}

/* Absolute value */
inline dd_real abs(const dd_real &a) {
  return (a.x[0] < 0.0) ? -a : a;
}

inline dd_real fabs(const dd_real &a) {
  return abs(a);
}

/* Round to Nearest integer */
inline dd_real nint(const dd_real &a) {
  double hi = qd::nint(a.x[0]);
  double lo;

  if (hi == a.x[0]) {
    /* High word is an integer already.  Round the low word.*/
    lo = qd::nint(a.x[1]);
    
    /* Renormalize. This is needed if x[0] = some integer, x[1] = 1/2.*/
    hi = qd::quick_two_sum(hi, lo, lo);
  } else {
    /* High word is not an integer. */
    lo = 0.0;
    if (std::abs(hi-a.x[0]) == 0.5 && a.x[1] < 0.0) {
      /* There is a tie in the high word, consult the low word 
         to break the tie. */
      hi -= 1.0;      /* NOTE: This does not cause INEXACT. */
    }
  }

  return dd_real(hi, lo);
}

inline dd_real floor(const dd_real &a) {
  double hi = std::floor(a.x[0]);
  double lo = 0.0;

  if (hi == a.x[0]) {
    /* High word is integer already.  Round the low word. */
    lo = std::floor(a.x[1]);
    hi = qd::quick_two_sum(hi, lo, lo);
  }

  return dd_real(hi, lo);
}

inline dd_real ceil(const dd_real &a) {
  double hi = std::ceil(a.x[0]);
  double lo = 0.0;

  if (hi == a.x[0]) {
    /* High word is integer already.  Round the low word. */
    lo = std::ceil(a.x[1]);
    hi = qd::quick_two_sum(hi, lo, lo);
  }

  return dd_real(hi, lo);
}

inline dd_real aint(const dd_real &a) {
  return (a.x[0] >= 0.0) ? floor(a) : ceil(a);
}

/* Cast to double. */
inline double to_double(const dd_real &a) {
  return a.x[0];
}

/* Cast to int. */
inline int to_int(const dd_real &a) {
  return static_cast<int>(a.x[0]);
}

/* Random number generator */
inline dd_real dd_real::rand() {
  return ddrand();
}

#endif /* _QD_DD_INLINE_H */
