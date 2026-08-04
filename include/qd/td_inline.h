/*
 * include/td_inline.h
 *
 * Contains small inline functions for the triple-double precision package.
 * Modeled on qd_inline.h.
 */
#ifndef _QD_TD_INLINE_H
#define _QD_TD_INLINE_H

#include <cmath>
#include <qd/inline.h>
#include <qd/qd_inline.h>   /* for qd::three_sum, qd::quick_three_accum */

#ifndef QD_INLINE
#define inline
#endif

/********** Constructors **********/
inline td_real::td_real(double x0, double x1, double x2) {
  x[0] = x0;
  x[1] = x1;
  x[2] = x2;
}

inline td_real::td_real(const double *xx) {
  x[0] = xx[0];
  x[1] = xx[1];
  x[2] = xx[2];
}

inline td_real::td_real(double x0) {
  x[0] = x0;
  x[1] = x[2] = 0.0;
}

inline td_real::td_real() {
  x[0] = 0.0;
  x[1] = 0.0;
  x[2] = 0.0;
}

inline td_real::td_real(const dd_real &a) {
  x[0] = a._hi();
  x[1] = a._lo();
  x[2] = 0.0;
}

inline td_real::td_real(const qd_real &a) {
  /* Renormalize the first three components after dropping a[3]. */
  double c0 = a[0], c1 = a[1], c2 = a[2], c3 = a[3];
  double s0, s1, s2, t;

  if (QD_ISINF(c0) || QD_ISNAN(c0)) {
    x[0] = c0; x[1] = c1; x[2] = c2;
    return;
  }

  /* Fold c3 upward. */
  s2 = qd::quick_two_sum(c2, c3, t);
  s1 = qd::quick_two_sum(c1, s2, s2);
  s0 = qd::quick_two_sum(c0, s1, s1);
  /* Now (s0, s1, s2) is renormalized; t is the discarded tail. */
  x[0] = s0;
  x[1] = s1;
  x[2] = s2 + t;
}

inline td_real::td_real(int i) {
  x[0] = static_cast<double>(i);
  x[1] = x[2] = 0.0;
}

/********** Accessors **********/
inline double td_real::operator[](int i) const {
  return x[i];
}

inline double &td_real::operator[](int i) {
  return x[i];
}

inline bool td_real::isnan() const {
  return QD_ISNAN(x[0]) || QD_ISNAN(x[1]) || QD_ISNAN(x[2]);
}

/********** Renormalization **********/
namespace qd {

/* Renormalize three components in place. */
inline void renorm3(double &c0, double &c1, double &c2) {
  double s0, s1, s2 = 0.0;

  if (QD_ISINF(c0)) return;

  s0 = qd::quick_two_sum(c1, c2, c2);
  c0 = qd::quick_two_sum(c0, s0, c1);

  s0 = c0;
  s1 = c1;
  if (s1 != 0.0) {
    s1 = qd::quick_two_sum(s1, c2, s2);
  } else {
    s0 = qd::quick_two_sum(s0, c2, s1);
  }

  c0 = s0;
  c1 = s1;
  c2 = s2;
}

/* Renormalize four components into three. */
inline void renorm3(double &c0, double &c1, double &c2, double &c3) {
  double s0, s1, s2 = 0.0;

  if (QD_ISINF(c0)) return;

  s0 = qd::quick_two_sum(c2, c3, c3);
  s0 = qd::quick_two_sum(c1, s0, c2);
  c0 = qd::quick_two_sum(c0, s0, c1);

  s0 = c0;
  s1 = c1;
  if (s1 != 0.0) {
    s1 = qd::quick_two_sum(s1, c2, s2);
    if (s2 != 0.0)
      s2 += c3;
    else
      s1 = qd::quick_two_sum(s1, c3, s2);
  } else {
    s0 = qd::quick_two_sum(s0, c2, s1);
    if (s1 != 0.0)
      s1 = qd::quick_two_sum(s1, c3, s2);
    else
      s0 = qd::quick_two_sum(s0, c3, s1);
  }

  c0 = s0;
  c1 = s1;
  c2 = s2;
}

}

inline void td_real::renorm() {
  qd::renorm3(x[0], x[1], x[2]);
}

inline void td_real::renorm(double &e) {
  qd::renorm3(x[0], x[1], x[2], e);
}

/********** Additions ************/

/* triple-double + double */
inline td_real operator+(const td_real &a, double b) {
  double c0, c1, c2;
  double e;

  c0 = qd::two_sum(a[0], b, e);
  c1 = qd::two_sum(a[1], e, e);
  c2 = qd::two_sum(a[2], e, e);

  qd::renorm3(c0, c1, c2, e);

  return td_real(c0, c1, c2);
}

/* triple-double + double-double */
inline td_real operator+(const td_real &a, const dd_real &b) {
  double s0, s1, s2;
  double t0, t1;

  s0 = qd::two_sum(a[0], b._hi(), t0);
  s1 = qd::two_sum(a[1], b._lo(), t1);

  s1 = qd::two_sum(s1, t0, t0);

  s2 = a[2];
  qd::three_sum(s2, t0, t1);

  /* Now (s0, s1, s2) plus tail (t0 + t1) */
  double tail = t0 + t1;
  qd::renorm3(s0, s1, s2, tail);
  return td_real(s0, s1, s2);
}

/* double + triple-double */
inline td_real operator+(double a, const td_real &b) {
  return (b + a);
}

/* double-double + triple-double */
inline td_real operator+(const dd_real &a, const td_real &b) {
  return (b + a);
}

inline td_real td_real::ieee_add(const td_real &a, const td_real &b) {
  int i, j, k;
  double s, t;
  double u, v;
  double xx[3] = {0.0, 0.0, 0.0};

  i = j = k = 0;
  if (std::abs(a[i]) > std::abs(b[j]))
    u = a[i++];
  else
    u = b[j++];
  if (std::abs(a[i]) > std::abs(b[j]))
    v = a[i++];
  else
    v = b[j++];

  u = qd::quick_two_sum(u, v, v);

  while (k < 3) {
    if (i >= 3 && j >= 3) {
      xx[k] = u;
      if (k < 2)
        xx[++k] = v;
      break;
    }

    if (i >= 3)
      t = b[j++];
    else if (j >= 3)
      t = a[i++];
    else if (std::abs(a[i]) > std::abs(b[j])) {
      t = a[i++];
    } else
      t = b[j++];

    s = qd::quick_three_accum(u, v, t);

    if (s != 0.0) {
      xx[k++] = s;
    }
  }

  /* add the rest. */
  for (k = i; k < 3; k++)
    xx[2] += a[k];
  for (k = j; k < 3; k++)
    xx[2] += b[k];

  qd::renorm3(xx[0], xx[1], xx[2]);
  return td_real(xx[0], xx[1], xx[2]);
}

inline td_real td_real::sloppy_add(const td_real &a, const td_real &b) {
  double s0, s1, s2;
  double t0, t1, t2;

  double v0, v1, v2;
  double u0, u1, u2;
  double w0, w1, w2;

  s0 = a[0] + b[0];
  s1 = a[1] + b[1];
  s2 = a[2] + b[2];

  v0 = s0 - a[0];
  v1 = s1 - a[1];
  v2 = s2 - a[2];

  u0 = s0 - v0;
  u1 = s1 - v1;
  u2 = s2 - v2;

  w0 = a[0] - u0;
  w1 = a[1] - u1;
  w2 = a[2] - u2;

  u0 = b[0] - v0;
  u1 = b[1] - v1;
  u2 = b[2] - v2;

  t0 = w0 + u0;
  t1 = w1 + u1;
  t2 = w2 + u2;

  s1 = qd::two_sum(s1, t0, t0);
  qd::three_sum(s2, t0, t1);
  t0 = t0 + t1 + t2;

  qd::renorm3(s0, s1, s2, t0);
  return td_real(s0, s1, s2);
}

// 2025-12-24(Wed) T.Kouya
// Branch free algorithm
//void Add3(const double x[3], const double y[3], double z[3]) {
//static inline void c_td_add_bf(double *a, double *b, double *c)
inline td_real td_real::bf_add(const td_real &a, const td_real &b) {
  double a0, b0, c0, d0, e0, f0;
  double a1, b1, c1, d1, e1, f1;
  double a2, b2, c2, d2, e2;
  double a3, b3, c3, d3;
  double c4;
  double c5, d5;
  double b6, c6;
  double a7, b7, c7;
  double b8, c8;

  a0 = a[0];
  b0 = b[0];
  c0 = a[1];
  d0 = b[1];
  e0 = a[2];
  f0 = b[2];
  a1 = qd::two_sum(a0, b0, b1);
  c1 = qd::two_sum(c0, d0, d1);
  e1 = qd::two_sum(e0, f0, f1);
  a2 = qd::quick_two_sum(a1, c1, c2);
  b2 = b1 + f1;
  d2 = qd::two_sum(d1, e1, e2);
  a3 = qd::quick_two_sum(a2, d2, d3);
  b3 = qd::two_sum(b2, c2, c3);
  c4 = c3 + e2;
  c5 = qd::two_sum(c4, d3, d5);
  b6 = qd::two_sum(b3, c5, c6);
  a7 = qd::quick_two_sum(a3, b6, b7);
  c7 = c6 + d5;
  b8 = qd::quick_two_sum(b7, c7, c8);

  return td_real(a7, b8, c8);
}

/* triple-double + triple-double */
inline td_real operator+(const td_real &a, const td_real &b) {
#ifndef QD_IEEE_ADD
  //return td_real::sloppy_add(a, b);
  return td_real::bf_add(a, b);
#else
  return td_real::ieee_add(a, b);
#endif
}

/********** Self-Additions ************/
inline td_real &td_real::operator+=(double a) {
  *this = *this + a;
  return *this;
}

inline td_real &td_real::operator+=(const dd_real &a) {
  *this = *this + a;
  return *this;
}

inline td_real &td_real::operator+=(const td_real &a) {
  *this = *this + a;
  return *this;
}

/********** Unary Minus **********/
inline td_real td_real::operator-() const {
  return td_real(-x[0], -x[1], -x[2]);
}

/********** Subtractions **********/
inline td_real operator-(const td_real &a, double b) {
  return (a + (-b));
}

inline td_real operator-(double a, const td_real &b) {
  return (a + (-b));
}

inline td_real operator-(const td_real &a, const dd_real &b) {
  return (a + (-b));
}

inline td_real operator-(const dd_real &a, const td_real &b) {
  return (a + (-b));
}

inline td_real operator-(const td_real &a, const td_real &b) {
  return (a + (-b));
}

/********** Self-Subtractions **********/
inline td_real &td_real::operator-=(double a) {
  return ((*this) += (-a));
}

inline td_real &td_real::operator-=(const dd_real &a) {
  return ((*this) += (-a));
}

inline td_real &td_real::operator-=(const td_real &a) {
  return ((*this) += (-a));
}


inline td_real operator*(double a, const td_real &b) {
  return (b * a);
}

inline td_real operator*(const dd_real &a, const td_real &b) {
  return (b * a);
}

inline td_real mul_pwr2(const td_real &a, double b) {
  return td_real(a[0] * b, a[1] * b, a[2] * b);
}

/********** Multiplications **********/
inline td_real operator*(const td_real &a, double b) {
  double p0, p1, p2;
  double q0, q1, q2;

  p0 = qd::two_prod(a[0], b, q0);   /* level 0/1 */
  p1 = qd::two_prod(a[1], b, q1);   /* level 1/2 */
  p2 = qd::two_prod(a[2], b, q2);   /* level 2/3 */

  /* Combine level-1 terms (q0, p1) and level-2 terms (p2, q1). */
  qd::three_sum(p1, p2, q0);

  /* After three_sum, p1 dominates level 1, p2 level 2, q0 smaller.
     Fold remaining smaller terms (q0, q1, q2) into a single tail. */
  double tail = q0 + q1 + q2;

  qd::renorm3(p0, p1, p2, tail);
  return td_real(p0, p1, p2);
}

/* triple-double * double-double */
/* a0 * b0                    0
        a0 * b1               1
        a1 * b0               2
             a1 * b1          3
             a2 * b0          4
                  a2 * b1     5  */
inline td_real operator*(const td_real &a, const dd_real &b) {
  double p0, p1, p2;
  double q0, q1, q2;

  p0 = qd::two_prod(a[0], b._hi(), q0);
  p1 = qd::two_prod(a[0], b._lo(), q1);
  p2 = qd::two_prod(a[1], b._hi(), q2);

  qd::three_sum(p1, p2, q0);

  /* q0 (now small from three_sum) plus all O(eps^2) and finer terms. */
  double tail = q0 + q1 + q2
              + a[1] * b._lo() + a[2] * b._hi() + a[2] * b._lo();

  qd::renorm3(p0, p1, p2, tail);
  return td_real(p0, p1, p2);
}

/* triple-double * triple-double (sloppy) */
/* a0 * b0                    0
        a0 * b1               1
        a1 * b0               2
             a0 * b2          3
             a1 * b1          4
             a2 * b0          5  */
inline td_real td_real::sloppy_mul(const td_real &a, const td_real &b) {
  double p0, p1, p2;
  double q0, q1, q2;

  p0 = qd::two_prod(a[0], b[0], q0);

  p1 = qd::two_prod(a[0], b[1], q1);
  p2 = qd::two_prod(a[1], b[0], q2);

  /* Accumulate (p1, p2, q0) at level eps. */
  qd::three_sum(p1, p2, q0);

  /* O(eps^2) and finer terms collected into a single double. */
  double tail = q1 + q2
              + a[0] * b[2] + a[1] * b[1] + a[2] * b[0];

  /* Now we have approximately p0 + p1 + p2 + (q0 + tail). */
  double s = q0 + tail;
  qd::renorm3(p0, p1, p2, s);
  return td_real(p0, p1, p2);
}

inline td_real td_real::accurate_mul(const td_real &a, const td_real &b) {
  double p0, p1, p2, p3, p4, p5;
  double q0, q1, q2, q3, q4, q5;
  double t0, t1;
  double s0, s1, s2;

  p0 = qd::two_prod(a[0], b[0], q0);

  p1 = qd::two_prod(a[0], b[1], q1);
  p2 = qd::two_prod(a[1], b[0], q2);

  p3 = qd::two_prod(a[0], b[2], q3);
  p4 = qd::two_prod(a[1], b[1], q4);
  p5 = qd::two_prod(a[2], b[0], q5);

  /* Accumulate level eps. */
  qd::three_sum(p1, p2, q0);

  /* Six-Three Sum  of p2, q1, q2, p3, p4, p5. */
  qd::three_sum(p2, q1, q2);
  qd::three_sum(p3, p4, p5);
  s0 = qd::two_sum(p2, p3, t0);
  s1 = qd::two_sum(q1, p4, t1);
  s2 = q2 + p5;
  s1 = qd::two_sum(s1, t0, t0);
  s2 += (t0 + t1);

  /* Higher-order tail (eps^3). */
  s1 += a[1] * b[2] + a[2] * b[1] + a[2] * b[2] + q3 + q4 + q5;

  /* Result: p0, p1, s0, with extra tail s1+s2 to be folded in. */
  double tail = s1 + s2;
  qd::renorm3(p0, p1, s0, tail);
  return td_real(p0, p1, s0);
}

// 2025-12-24(Wed) T.Kouya
// Branch free algorithm
// void Mul3(const double x[3], const double y[3], double z[3]) {
//static inline void c_td_mul_bf(double *a, double *b, double *c)
inline td_real td_real::bf_mul(const td_real &a, const td_real &b) {
  double a0, b0, c0, d0, e0, f0, g0, h0, i0;
  double c1, d1, e1, f1, g1, h1, i1;
  double b2, c2, g2;
  double a3, b3, c3, d3, e3, g3;
  double c4, e4;
  double b5, c5;
  double a6, b6;
  double b7, c7;

  a0 = qd::two_prod(a[0], b[0], b0);
  c0 = qd::two_prod(a[0], b[1], e0);
  d0 = qd::two_prod(a[1], b[0], f0);
  g0 = a[0] * b[2];
  h0 = a[1] * b[1];
  i0 = a[2] * b[0];
  c1 = qd::two_sum(c0, d0, d1);
  e1 = qd::two_sum(e0, f0, f1);
  g1 = qd::two_sum(g0, i0, i1);
  b2 = qd::two_sum(b0, c1, c2);
  g2 = qd::two_sum(g1, h0, h1);
  a3 = qd::quick_two_sum(a0, b2, b3);
  c3 = qd::two_sum(c2, d1, d3);
  e3 = qd::two_sum(e1, g2, g3);
  c4 = qd::two_sum(c3, e3, e4);
  b5 = qd::quick_two_sum(b3, c4, c5);
  a6 = qd::quick_two_sum(a3, b5, b6);
  b7 = qd::quick_two_sum(b6, c5, c7);

  return td_real(a6, b7, c7);
}

inline td_real operator*(const td_real &a, const td_real &b) {
#ifdef QD_SLOPPY_MUL
  //return td_real::sloppy_mul(a, b);
  return td_real::bf_mul(a, b);
#else
  return td_real::accurate_mul(a, b);
  //return td_real::bf_mul(a, b);
#endif
}

// 2026-08-04 T.Kouya
// Branch free algorithm: triple-word FMA,  z = a * b + c.
//
// Every term of the exact product a*b and every word of c is assigned to
// the magnitude level it belongs to (level 0 ~ 1, level 1 ~ eps,
// level 2 ~ eps^2).  The levels are accumulated with two_sum, the errors
// spilling one level down, and the three level accumulators are
// renormalized by the same quick_two_sum network used by td_real::bf_mul.
// No branch anywhere, and a*b is never renormalized on its own.
inline td_real tw_fma(const td_real &a, const td_real &b, const td_real &c) {
  double a0, b0, c0, d0, e0, f0, g0;
  double a1, b1, c1, d1, e1, f1;
  double a2, b2, c2;
  double a3, b3;
  double a4, b4;
  double a5, b5;
  double a6, b6;

  /* --- products ------------------------------------------------- */
  a0 = qd::two_prod(a[0], b[0], b0);         /* level 0 / level 1 */
  c0 = qd::two_prod(a[0], b[1], e0);         /* level 1 / level 2 */
  d0 = qd::two_prod(a[1], b[0], f0);         /* level 1 / level 2 */
  g0 = a[0] * b[2] + a[1] * b[1] + a[2] * b[0];   /* level 2 */

  /* --- level 0 : leading product word + c[0] --------------------- */
  a1 = qd::two_sum(a0, c[0], b1);            /* b1 spills to level 1 */

  /* --- level 1 : b0, c0, d0, c[1], b1 ---------------------------- */
  c1 = qd::two_sum(c0, d0, d1);
  e1 = qd::two_sum(b0, c[1], f1);
  a2 = qd::two_sum(c1, e1, b2);
  a3 = qd::two_sum(a2, b1, b3);              /* level-1 accumulator a3 */

  /* --- level 2 : all spills and the eps^2 product terms ---------- */
  c2 = ((e0 + f0) + (g0 + c[2])) + ((d1 + f1) + (b2 + b3));

  /* --- branch-free renormalization (as in td_real::bf_mul) ------- */
  a4 = qd::quick_two_sum(a3, c2, b4);        /* level 1 + level 2 */
  a5 = qd::quick_two_sum(a1, a4, b5);        /* level 0 + level 1 */
  a6 = qd::quick_two_sum(b5, b4, b6);

  return td_real(a5, a6, b6);
}

/* triple-word FMA with a plain double multiplier:  a * b + c. */
inline td_real tw_fma(const td_real &a, double b, const td_real &c) {
  double a0, b0, c0, d0, e0;
  double a1, b1, c1, d1;
  double a2, b2;
  double a3, b3, c3;
  double a4, b4;
  double a5, b5;
  double a6, b6;

  a0 = qd::two_prod(a[0], b, b0);            /* level 0 / level 1 */
  c0 = qd::two_prod(a[1], b, d0);            /* level 1 / level 2 */
  e0 = a[2] * b;                             /* level 2 */

  a1 = qd::two_sum(a0, c[0], b1);            /* level 0, spill b1 */

  c1 = qd::two_sum(b0, c0, d1);
  a2 = qd::two_sum(c[1], b1, b2);
  a3 = qd::two_sum(c1, a2, b3);              /* level-1 accumulator a3 */

  c3 = ((d0 + e0) + c[2]) + ((d1 + b2) + b3);

  /* branch-free renormalization */
  a4 = qd::quick_two_sum(a3, c3, b4);        /* level 1 + level 2 */
  a5 = qd::quick_two_sum(a1, a4, b5);        /* level 0 + level 1 */
  a6 = qd::quick_two_sum(b5, b4, b6);

  return td_real(a5, a6, b6);
}

inline td_real fma(const td_real &a, const td_real &b, const td_real &c) {
  return tw_fma(a, b, c);
}

inline td_real fma(const td_real &a, double b, const td_real &c) {
  return tw_fma(a, b, c);
}

/* triple-double ^ 2  = (x0 + x1 + x2) ^ 2
                      = x0^2 + 2 x0 x1 + (2 x0 x2 + x1^2) + 2 x1 x2 + x2^2 */
inline td_real sqr(const td_real &a) {
  double p0, p1;
  double q0, q1;

  p0 = qd::two_sqr(a[0], q0);                    /* a0^2  level 0/1 */
  p1 = qd::two_prod(2.0 * a[0], a[1], q1);       /* 2 a0 a1  level 1/2 */

  /* Level 1 terms: p1 (eps), q0 (eps).  q1 is eps^2. */
  qd::three_sum(p1, q0, q1);

  /* Level 2 and finer terms */
  double t = 2.0 * a[0] * a[2] + a[1] * a[1]
           + 2.0 * a[1] * a[2] + a[2] * a[2];

  /* q1 (small from three_sum) plus level-2+ terms. */
  q1 += t;

  qd::renorm3(p0, p1, q0, q1);
  return td_real(p0, p1, q0);
}

/********** Self-Multiplication **********/
inline td_real &td_real::operator*=(double a) {
  *this = (*this * a);
  return *this;
}

inline td_real &td_real::operator*=(const dd_real &a) {
  *this = (*this * a);
  return *this;
}

inline td_real &td_real::operator*=(const td_real &a) {
  *this = *this * a;
  return *this;
}

inline td_real operator/(const td_real &a, const dd_real &b) {
#ifdef QD_SLOPPY_DIV
  return td_real::sloppy_div(a, b);
#else
  return td_real::accurate_div(a, b);
#endif
}

inline td_real operator/(const td_real &a, const td_real &b) {
#ifdef QD_SLOPPY_DIV
  return td_real::sloppy_div(a, b);
#elif defined(QD_NO_FMA_DIV)
  return td_real::accurate_div(a, b);
#else
  return td_real::fma_div(a, b);
#endif
}

/* double / triple-double */
inline td_real operator/(double a, const td_real &b) {
  return td_real(a) / b;
}

/* double-double / triple-double */
inline td_real operator/(const dd_real &a, const td_real &b) {
  return td_real(a) / b;
}

/********** Self-Divisions **********/
inline td_real &td_real::operator/=(double a) {
  *this = (*this / a);
  return *this;
}

inline td_real &td_real::operator/=(const dd_real &a) {
  *this = (*this / a);
  return *this;
}

inline td_real &td_real::operator/=(const td_real &a) {
  *this = (*this / a);
  return *this;
}


/********** Exponentiation **********/
inline td_real td_real::operator^(int n) const {
  return pow(*this, n);
}

/********** Miscellaneous **********/
inline td_real abs(const td_real &a) {
  return (a[0] < 0.0) ? -a : a;
}

inline td_real fabs(const td_real &a) {
  return abs(a);
}

/* Quick version.  May be off by one when td is very close
   to the middle of two integers.                         */
inline td_real quick_nint(const td_real &a) {
  td_real r = td_real(qd::nint(a[0]), qd::nint(a[1]), qd::nint(a[2]));
  r.renorm();
  return r;
}

/*********** Assignments ************/
inline td_real &td_real::operator=(double a) {
  x[0] = a;
  x[1] = x[2] = 0.0;
  return *this;
}

inline td_real &td_real::operator=(const dd_real &a) {
  x[0] = a._hi();
  x[1] = a._lo();
  x[2] = 0.0;
  return *this;
}

inline td_real &td_real::operator=(const qd_real &a) {
  *this = td_real(a);
  return *this;
}

/********** Equality Comparison **********/
inline bool operator==(const td_real &a, double b) {
  return (a[0] == b && a[1] == 0.0 && a[2] == 0.0);
}

inline bool operator==(double a, const td_real &b) {
  return (b == a);
}

inline bool operator==(const td_real &a, const dd_real &b) {
  return (a[0] == b._hi() && a[1] == b._lo() && a[2] == 0.0);
}

inline bool operator==(const dd_real &a, const td_real &b) {
  return (b == a);
}

inline bool operator==(const td_real &a, const td_real &b) {
  return (a[0] == b[0] && a[1] == b[1] && a[2] == b[2]);
}


/********** Less-Than Comparison ***********/
inline bool operator<(const td_real &a, double b) {
  return (a[0] < b || (a[0] == b && a[1] < 0.0));
}

inline bool operator<(double a, const td_real &b) {
  return (b > a);
}

inline bool operator<(const td_real &a, const dd_real &b) {
  return (a[0] < b._hi() ||
          (a[0] == b._hi() && (a[1] < b._lo() ||
                            (a[1] == b._lo() && a[2] < 0.0))));
}

inline bool operator<(const dd_real &a, const td_real &b) {
  return (b > a);
}

inline bool operator<(const td_real &a, const td_real &b) {
  return (a[0] < b[0] ||
          (a[0] == b[0] && (a[1] < b[1] ||
                            (a[1] == b[1] && a[2] < b[2]))));
}

/********** Greater-Than Comparison ***********/
inline bool operator>(const td_real &a, double b) {
  return (a[0] > b || (a[0] == b && a[1] > 0.0));
}

inline bool operator>(double a, const td_real &b) {
  return (b < a);
}

inline bool operator>(const td_real &a, const dd_real &b) {
  return (a[0] > b._hi() ||
          (a[0] == b._hi() && (a[1] > b._lo() ||
                            (a[1] == b._lo() && a[2] > 0.0))));
}

inline bool operator>(const dd_real &a, const td_real &b) {
  return (b < a);
}

inline bool operator>(const td_real &a, const td_real &b) {
  return (a[0] > b[0] ||
          (a[0] == b[0] && (a[1] > b[1] ||
                            (a[1] == b[1] && a[2] > b[2]))));
}


/********** Less-Than-Or-Equal-To Comparison **********/
inline bool operator<=(const td_real &a, double b) {
  return (a[0] < b || (a[0] == b && a[1] <= 0.0));
}

inline bool operator<=(double a, const td_real &b) {
  return (b >= a);
}

inline bool operator<=(const td_real &a, const dd_real &b) {
  return (a[0] < b._hi() ||
          (a[0] == b._hi() && (a[1] < b._lo() ||
                            (a[1] == b._lo() && a[2] <= 0.0))));
}

inline bool operator<=(const dd_real &a, const td_real &b) {
  return (b >= a);
}

inline bool operator<=(const td_real &a, const td_real &b) {
  return (a[0] < b[0] ||
          (a[0] == b[0] && (a[1] < b[1] ||
                            (a[1] == b[1] && a[2] <= b[2]))));
}

/********** Greater-Than-Or-Equal-To Comparison **********/
inline bool operator>=(const td_real &a, double b) {
  return (a[0] > b || (a[0] == b && a[1] >= 0.0));
}

inline bool operator>=(double a, const td_real &b) {
  return (b <= a);
}

inline bool operator>=(const td_real &a, const dd_real &b) {
  return (a[0] > b._hi() ||
          (a[0] == b._hi() && (a[1] > b._lo() ||
                            (a[1] == b._lo() && a[2] >= 0.0))));
}

inline bool operator>=(const dd_real &a, const td_real &b) {
  return (b <= a);
}

inline bool operator>=(const td_real &a, const td_real &b) {
  return (a[0] > b[0] ||
          (a[0] == b[0] && (a[1] > b[1] ||
                            (a[1] == b[1] && a[2] >= b[2]))));
}

/********** Not-Equal-To Comparison **********/
inline bool operator!=(const td_real &a, double b) {
  return !(a == b);
}

inline bool operator!=(double a, const td_real &b) {
  return !(a == b);
}

inline bool operator!=(const td_real &a, const dd_real &b) {
  return !(a == b);
}

inline bool operator!=(const dd_real &a, const td_real &b) {
  return !(a == b);
}

inline bool operator!=(const td_real &a, const td_real &b) {
  return !(a == b);
}



inline td_real aint(const td_real &a) {
  return (a[0] >= 0) ? floor(a) : ceil(a);
}

inline bool td_real::is_zero() const {
  return (x[0] == 0.0);
}

inline bool td_real::is_one() const {
  return (x[0] == 1.0 && x[1] == 0.0 && x[2] == 0.0);
}

inline bool td_real::is_positive() const {
  return (x[0] > 0.0);
}

inline bool td_real::is_negative() const {
  return (x[0] < 0.0);
}

inline dd_real to_dd_real(const td_real &a) {
  return dd_real(a[0], a[1]);
}

inline qd_real to_qd_real(const td_real &a) {
  return qd_real(a[0], a[1], a[2], 0.0);
}

inline td_real to_td_real(const qd_real &a) {
  return td_real(a);
}

inline double to_double(const td_real &a) {
  return a[0];
}

inline int to_int(const td_real &a) {
  return static_cast<int>(a[0]);
}

inline td_real inv(const td_real &td) {
  return 1.0 / td;
}

inline td_real max(const td_real &a, const td_real &b) {
  return (a > b) ? a : b;
}

inline td_real max(const td_real &a, const td_real &b,
                   const td_real &c) {
  return (a > b) ? ((a > c) ? a : c) : ((b > c) ? b : c);
}

inline td_real min(const td_real &a, const td_real &b) {
  return (a < b) ? a : b;
}

inline td_real min(const td_real &a, const td_real &b,
                   const td_real &c) {
  return (a < b) ? ((a < c) ? a : c) : ((b < c) ? b : c);
}

/* Random number generator */
inline td_real td_real::rand() {
  return tdrand();
}

inline td_real ldexp(const td_real &a, int n) {
  return td_real(std::ldexp(a[0], n), std::ldexp(a[1], n),
                 std::ldexp(a[2], n));
}

#endif /* _QD_TD_INLINE_H */
