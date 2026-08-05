/*
 * include/qd_inline.h
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2000-2001
 *
 * Contains small functions (suitable for inlining) in the quad-double
 * arithmetic package.
 */
#ifndef _QD_QD_INLINE_H
#define _QD_QD_INLINE_H

#include <cmath>
#include <qd/inline.h>

#ifndef QD_INLINE
#define inline
#endif

/********** Constructors **********/
inline qd_real::qd_real(double x0, double x1, double x2, double x3) {
  x[0] = x0;
  x[1] = x1;
  x[2] = x2;
  x[3] = x3;
}

inline qd_real::qd_real(const double *xx) {
  x[0] = xx[0];
  x[1] = xx[1];
  x[2] = xx[2];
  x[3] = xx[3];
}

inline qd_real::qd_real(double x0) {
  x[0] = x0;
  x[1] = x[2] = x[3] = 0.0;
}

inline qd_real::qd_real() {
	x[0] = 0.0; 
	x[1] = 0.0; 
	x[2] = 0.0; 
	x[3] = 0.0; 
}

inline qd_real::qd_real(const dd_real &a) {
  x[0] = a._hi();
  x[1] = a._lo();
  x[2] = x[3] = 0.0;
}

inline qd_real::qd_real(int i) {
  x[0] = static_cast<double>(i);
  x[1] = x[2] = x[3] = 0.0;
}

/********** Accessors **********/
inline double qd_real::operator[](int i) const {
  return x[i];
}

inline double &qd_real::operator[](int i) {
  return x[i];
}

inline bool qd_real::isnan() const {
  return QD_ISNAN(x[0]) || QD_ISNAN(x[1]) || QD_ISNAN(x[2]) || QD_ISNAN(x[3]);
}

/********** Renormalization **********/
namespace qd {
inline void quick_renorm(double &c0, double &c1, 
                         double &c2, double &c3, double &c4) {
  double t0, t1, t2, t3;
  double s;
  s  = qd::quick_two_sum(c3, c4, t3);
  s  = qd::quick_two_sum(c2, s , t2);
  s  = qd::quick_two_sum(c1, s , t1);
  c0 = qd::quick_two_sum(c0, s , t0);

  s  = qd::quick_two_sum(t2, t3, t2);
  s  = qd::quick_two_sum(t1, s , t1);
  c1 = qd::quick_two_sum(t0, s , t0);

  s  = qd::quick_two_sum(t1, t2, t1);
  c2 = qd::quick_two_sum(t0, s , t0);
  
  c3 = t0 + t1;
}

inline void renorm(double &c0, double &c1, 
                   double &c2, double &c3) {
  double s0, s1, s2 = 0.0, s3 = 0.0;

  if (QD_ISINF(c0)) return;

  s0 = qd::quick_two_sum(c2, c3, c3);
  s0 = qd::quick_two_sum(c1, s0, c2);
  c0 = qd::quick_two_sum(c0, s0, c1);

  s0 = c0;
  s1 = c1;
  if (s1 != 0.0) {
    s1 = qd::quick_two_sum(s1, c2, s2);
    if (s2 != 0.0)
      s2 = qd::quick_two_sum(s2, c3, s3);
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
  c3 = s3;
}

inline void renorm(double &c0, double &c1, 
                   double &c2, double &c3, double &c4) {
  double s0, s1, s2 = 0.0, s3 = 0.0;

  if (QD_ISINF(c0)) return;

  s0 = qd::quick_two_sum(c3, c4, c4);
  s0 = qd::quick_two_sum(c2, s0, c3);
  s0 = qd::quick_two_sum(c1, s0, c2);
  c0 = qd::quick_two_sum(c0, s0, c1);

  s0 = c0;
  s1 = c1;

  s0 = qd::quick_two_sum(c0, c1, s1);
  if (s1 != 0.0) {
    s1 = qd::quick_two_sum(s1, c2, s2);
    if (s2 != 0.0) {
      s2 = qd::quick_two_sum(s2, c3, s3);
      if (s3 != 0.0)
        s3 += c4;
      else
        s2 += c4;
    } else {
      s1 = qd::quick_two_sum(s1, c3, s2);
      if (s2 != 0.0)
        s2 = qd::quick_two_sum(s2, c4, s3);
      else
        s1 = qd::quick_two_sum(s1, c4, s2);
    }
  } else {
    s0 = qd::quick_two_sum(s0, c2, s1);
    if (s1 != 0.0) {
      s1 = qd::quick_two_sum(s1, c3, s2);
      if (s2 != 0.0)
        s2 = qd::quick_two_sum(s2, c4, s3);
      else
        s1 = qd::quick_two_sum(s1, c4, s2);
    } else {
      s0 = qd::quick_two_sum(s0, c3, s1);
      if (s1 != 0.0)
        s1 = qd::quick_two_sum(s1, c4, s2);
      else
        s0 = qd::quick_two_sum(s0, c4, s1);
    }
  }

  c0 = s0;
  c1 = s1;
  c2 = s2;
  c3 = s3;
}
}

inline void qd_real::renorm() {
  qd::renorm(x[0], x[1], x[2], x[3]);
}

inline void qd_real::renorm(double &e) {
  qd::renorm(x[0], x[1], x[2], x[3], e);
}


/********** Additions ************/
namespace qd {

inline void three_sum(double &a, double &b, double &c) {
  double t1, t2, t3;
  t1 = qd::two_sum(a, b, t2);
  a  = qd::two_sum(c, t1, t3);
  b  = qd::two_sum(t2, t3, c);
}

inline void three_sum2(double &a, double &b, double &c) {
  double t1, t2, t3;
  t1 = qd::two_sum(a, b, t2);
  a  = qd::two_sum(c, t1, t3);
  b = t2 + t3;
}

}

/* quad-double + double */
inline qd_real operator+(const qd_real &a, double b) {
  double c0, c1, c2, c3;
  double e;

  c0 = qd::two_sum(a[0], b, e);
  c1 = qd::two_sum(a[1], e, e);
  c2 = qd::two_sum(a[2], e, e);
  c3 = qd::two_sum(a[3], e, e);

  qd::renorm(c0, c1, c2, c3, e);

  return qd_real(c0, c1, c2, c3);
}

/* quad-double + double-double */
inline qd_real operator+(const qd_real &a, const dd_real &b) {

  double s0, s1, s2, s3;
  double t0, t1;

  s0 = qd::two_sum(a[0], b._hi(), t0);
  s1 = qd::two_sum(a[1], b._lo(), t1);

  s1 = qd::two_sum(s1, t0, t0);

  s2 = a[2];
  qd::three_sum(s2, t0, t1);

  s3 = qd::two_sum(t0, a[3], t0);
  t0 += t1;

  qd::renorm(s0, s1, s2, s3, t0);
  return qd_real(s0, s1, s2, s3);
}


/* double + quad-double */
inline qd_real operator+(double a, const qd_real &b) {
  return (b + a);
}

/* double-double + quad-double */
inline qd_real operator+(const dd_real &a, const qd_real &b) {
  return (b + a);
}

namespace qd {

/* s = quick_three_accum(a, b, c) adds c to the dd-pair (a, b).
 * If the result does not fit in two doubles, then the sum is 
 * output into s and (a,b) contains the remainder.  Otherwise
 * s is zero and (a,b) contains the sum. */
inline double quick_three_accum(double &a, double &b, double c) {
  double s;
  bool za, zb;

  s = qd::two_sum(b, c, b);
  s = qd::two_sum(a, s, a);

  za = (a != 0.0);
  zb = (b != 0.0);

  if (za && zb)
    return s;

  if (!zb) {
    b = a;
    a = s;
  } else {
    a = s;
  }

  return 0.0;
}

}

inline qd_real qd_real::ieee_add(const qd_real &a, const qd_real &b) {
  int i, j, k;
  double s, t;
  double u, v;   /* double-length accumulator */
  double x[4] = {0.0, 0.0, 0.0, 0.0};
  
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
  
  while (k < 4) {
    if (i >= 4 && j >= 4) {
      x[k] = u;
      if (k < 3)
        x[++k] = v;
      break;
    }

    if (i >= 4)
      t = b[j++];
    else if (j >= 4)
      t = a[i++];
    else if (std::abs(a[i]) > std::abs(b[j])) {
      t = a[i++];
    } else
      t = b[j++];

    s = qd::quick_three_accum(u, v, t);

    if (s != 0.0) {
      x[k++] = s;
    }
  }

  /* add the rest. */
  for (k = i; k < 4; k++)
    x[3] += a[k];
  for (k = j; k < 4; k++)
    x[3] += b[k];

  qd::renorm(x[0], x[1], x[2], x[3]);
  return qd_real(x[0], x[1], x[2], x[3]);
}

inline qd_real qd_real::sloppy_add(const qd_real &a, const qd_real &b) {
  /*
  double s0, s1, s2, s3;
  double t0, t1, t2, t3;
  
  s0 = qd::two_sum(a[0], b[0], t0);
  s1 = qd::two_sum(a[1], b[1], t1);
  s2 = qd::two_sum(a[2], b[2], t2);
  s3 = qd::two_sum(a[3], b[3], t3);

  s1 = qd::two_sum(s1, t0, t0);
  qd::three_sum(s2, t0, t1);
  qd::three_sum2(s3, t0, t2);
  t0 = t0 + t1 + t3;

  qd::renorm(s0, s1, s2, s3, t0);
  return qd_real(s0, s1, s2, s3, t0);
  */

  /* Same as above, but addition re-organized to minimize
     data dependency ... unfortunately some compilers are
     not very smart to do this automatically */
  double s0, s1, s2, s3;
  double t0, t1, t2, t3;

  double v0, v1, v2, v3;
  double u0, u1, u2, u3;
  double w0, w1, w2, w3;

  s0 = a[0] + b[0];
  s1 = a[1] + b[1];
  s2 = a[2] + b[2];
  s3 = a[3] + b[3];

  v0 = s0 - a[0];
  v1 = s1 - a[1];
  v2 = s2 - a[2];
  v3 = s3 - a[3];

  u0 = s0 - v0;
  u1 = s1 - v1;
  u2 = s2 - v2;
  u3 = s3 - v3;

  w0 = a[0] - u0;
  w1 = a[1] - u1;
  w2 = a[2] - u2;
  w3 = a[3] - u3;

  u0 = b[0] - v0;
  u1 = b[1] - v1;
  u2 = b[2] - v2;
  u3 = b[3] - v3;

  t0 = w0 + u0;
  t1 = w1 + u1;
  t2 = w2 + u2;
  t3 = w3 + u3;

  s1 = qd::two_sum(s1, t0, t0);
  qd::three_sum(s2, t0, t1);
  qd::three_sum2(s3, t0, t2);
  t0 = t0 + t1 + t3;

  /* renormalize */
  qd::renorm(s0, s1, s2, s3, t0);
  return qd_real(s0, s1, s2, s3);
}

// 2025-12-25(Wed) T.Kouya
// Branch free algorithm
// void Add4(const double x[4], const double y[4], double z[4]) {
//static inline void c_qd_add_bf(const double *a, const double *b, double *c)
inline qd_real qd_real::bf_add(const qd_real &a, const qd_real &b) {
	double a0 , b0 , c0 , d0 , e0, f0, g0, h0;
	double a1 , b1 , c1 , d1 , e1, f1, g1, h1;
	double a2 , b2 , c2 , d2 , e2, f2, g2;
	double a3 , b3 , c3 , d3 , e3, f3, g3;
	double a4 , b4 , c4 , d4 , e4;
	double a5 , b5 , c5 , d5 , e5;
	double a6 , b6 , c6 , d6 , e6;
	double a7 , b7 , c7 , d7 , e7;
	double a8 , b8 , c8 , d8 , e8;
	double a9 , b9 , c9 , d9 ;
	double a10, b10, c10, d10;
	double a11, b11, c11, d11;
	double a12, b12, c12, d12;

	a0 = a[0];
  b0 = b[0];
  c0 = a[1];
  d0 = b[1];
  e0 = a[2];
  f0 = b[2];
  g0 = a[3];
  h0 = b[3];
  a1 = qd::two_sum(a0, b0, b1);
  c1 = qd::two_sum(c0, d0, d1);
  e1 = qd::two_sum(e0, f0, f1);
  g1 = qd::two_sum(g0, h0, h1);
  a2 = qd::quick_two_sum(a1, c1, c2);
  b2 = b1 + h1;
  d2 = qd::two_sum(d1, e1, e2);
  f2 = qd::two_sum(f1, g1, g2);
  b3 = qd::two_sum(b2, g2, g3);
  c3 = qd::quick_two_sum(c2, d2, d3);
  e3 = qd::two_sum(e2, f2, f3);
  a4 = qd::quick_two_sum(a2, c3, c4);
  d4 = qd::quick_two_sum(d3, e3, e4);
  b5 = qd::two_sum(b3, d4, d5);
  e5 = e4 + f3;
  b6 = qd::two_sum(b5, c4, c6);
  d6 = qd::two_sum(d5, e5, e6);
  a7 = qd::quick_two_sum(a4, b6, b7);
  c7 = qd::quick_two_sum(c6, d6, d7);
  e8 = e6 + g3;
  b8 = qd::quick_two_sum(b7, c7, c8);
  d9 = d7 + e8;
  a10 = qd::quick_two_sum(a7, b8, b10);
  c10 = qd::quick_two_sum(c8, d9, d10);
  b11 = qd::quick_two_sum(b10, c10, c11);
  c12 = qd::quick_two_sum(c11, d10, d12);

	//c[0] = a10;
	//c[1] = b11;
	//c[2] = c12;
	//c[3] = d12;
  return qd_real(a10, b11, c12, d12);
}

/* quad-double + quad-double */
inline qd_real operator+(const qd_real &a, const qd_real &b) {
#ifndef QD_IEEE_ADD
  //return qd_real::sloppy_add(a, b);
  return qd_real::bf_add(a, b);
#else
  return qd_real::ieee_add(a, b);
#endif
}



/********** Self-Additions ************/
/* quad-double += double */
inline qd_real &qd_real::operator+=(double a) {
  *this = *this + a;
  return *this;
}

/* quad-double += double-double */
inline qd_real &qd_real::operator+=(const dd_real &a) {
  *this = *this + a;
  return *this;
}

/* quad-double += quad-double */
inline qd_real &qd_real::operator+=(const qd_real &a) {
  *this = *this + a;
  return *this;
}

/********** Unary Minus **********/
inline qd_real qd_real::operator-() const {
  return qd_real(-x[0], -x[1], -x[2], -x[3]);
}

/********** Subtractions **********/
inline qd_real operator-(const qd_real &a, double b) {
  return (a + (-b));
}

inline qd_real operator-(double a, const qd_real &b) {
  return (a + (-b));
}

inline qd_real operator-(const qd_real &a, const dd_real &b) {
  return (a + (-b));
}

inline qd_real operator-(const dd_real &a, const qd_real &b) {
  return (a + (-b));
}

inline qd_real operator-(const qd_real &a, const qd_real &b) {
  return (a + (-b));
}

/********** Self-Subtractions **********/
inline qd_real &qd_real::operator-=(double a) {
  return ((*this) += (-a));
}

inline qd_real &qd_real::operator-=(const dd_real &a) {
  return ((*this) += (-a));
}

inline qd_real &qd_real::operator-=(const qd_real &a) {
  return ((*this) += (-a));
}


inline qd_real operator*(double a, const qd_real &b) {
  return (b * a);
}

inline qd_real operator*(const dd_real &a, const qd_real &b) {
  return (b * a);
}

inline qd_real mul_pwr2(const qd_real &a, double b) {
  return qd_real(a[0] * b, a[1] * b, a[2] * b, a[3] * b);
}

/********** Multiplications **********/
inline qd_real operator*(const qd_real &a, double b) {
  double p0, p1, p2, p3;
  double q0, q1, q2;
  double s0, s1, s2, s3, s4;

  p0 = qd::two_prod(a[0], b, q0);
  p1 = qd::two_prod(a[1], b, q1);
  p2 = qd::two_prod(a[2], b, q2);
  p3 = a[3] * b;

  s0 = p0;

  s1 = qd::two_sum(q0, p1, s2);

  qd::three_sum(s2, q1, p2);

  qd::three_sum2(q1, q2, p3);
  s3 = q1;

  s4 = q2 + p2;

  qd::renorm(s0, s1, s2, s3, s4);
  return qd_real(s0, s1, s2, s3);

}

/* quad-double * double-double */
/* a0 * b0                        0
        a0 * b1                   1
        a1 * b0                   2
             a1 * b1              3
             a2 * b0              4
                  a2 * b1         5
                  a3 * b0         6
                       a3 * b1    7 */
inline qd_real operator*(const qd_real &a, const dd_real &b) {
  double p0, p1, p2, p3, p4;
  double q0, q1, q2, q3, q4;
  double s0, s1, s2;
  double t0, t1;

  p0 = qd::two_prod(a[0], b._hi(), q0);
  p1 = qd::two_prod(a[0], b._lo(), q1);
  p2 = qd::two_prod(a[1], b._hi(), q2);
  p3 = qd::two_prod(a[1], b._lo(), q3);
  p4 = qd::two_prod(a[2], b._hi(), q4);
  
  qd::three_sum(p1, p2, q0);
  
  /* Five-Three-Sum */
  qd::three_sum(p2, p3, p4);
  q1 = qd::two_sum(q1, q2, q2);
  s0 = qd::two_sum(p2, q1, t0);
  s1 = qd::two_sum(p3, q2, t1);
  s1 = qd::two_sum(s1, t0, t0);
  s2 = t0 + t1 + p4;
  p2 = s0;

  p3 = a[2] * b._hi() + a[3] * b._lo() + q3 + q4;
  qd::three_sum2(p3, q0, s1);
  p4 = q0 + s2;

  qd::renorm(p0, p1, p2, p3, p4);
  return qd_real(p0, p1, p2, p3);
}

/* quad-double * quad-double */
/* a0 * b0                    0
        a0 * b1               1
        a1 * b0               2
             a0 * b2          3
             a1 * b1          4
             a2 * b0          5
                  a0 * b3     6
                  a1 * b2     7
                  a2 * b1     8
                  a3 * b0     9  */
inline qd_real qd_real::sloppy_mul(const qd_real &a, const qd_real &b) {
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

  /* Start Accumulation */
  qd::three_sum(p1, p2, q0);

  /* Six-Three Sum  of p2, q1, q2, p3, p4, p5. */
  qd::three_sum(p2, q1, q2);
  qd::three_sum(p3, p4, p5);
  /* compute (s0, s1, s2) = (p2, q1, q2) + (p3, p4, p5). */
  s0 = qd::two_sum(p2, p3, t0);
  s1 = qd::two_sum(q1, p4, t1);
  s2 = q2 + p5;
  s1 = qd::two_sum(s1, t0, t0);
  s2 += (t0 + t1);

  /* O(eps^3) order terms */
  s1 += a[0]*b[3] + a[1]*b[2] + a[2]*b[1] + a[3]*b[0] + q0 + q3 + q4 + q5;
  qd::renorm(p0, p1, s0, s1, s2);
  return qd_real(p0, p1, s0, s1);
}

inline qd_real qd_real::accurate_mul(const qd_real &a, const qd_real &b) {
  double p0, p1, p2, p3, p4, p5;
  double q0, q1, q2, q3, q4, q5;
  double p6, p7, p8, p9;
  double q6, q7, q8, q9;
  double r0, r1;
  double t0, t1;
  double s0, s1, s2;

  p0 = qd::two_prod(a[0], b[0], q0);

  p1 = qd::two_prod(a[0], b[1], q1);
  p2 = qd::two_prod(a[1], b[0], q2);

  p3 = qd::two_prod(a[0], b[2], q3);
  p4 = qd::two_prod(a[1], b[1], q4);
  p5 = qd::two_prod(a[2], b[0], q5);

  /* Start Accumulation */
  qd::three_sum(p1, p2, q0);

  /* Six-Three Sum  of p2, q1, q2, p3, p4, p5. */
  qd::three_sum(p2, q1, q2);
  qd::three_sum(p3, p4, p5);
  /* compute (s0, s1, s2) = (p2, q1, q2) + (p3, p4, p5). */
  s0 = qd::two_sum(p2, p3, t0);
  s1 = qd::two_sum(q1, p4, t1);
  s2 = q2 + p5;
  s1 = qd::two_sum(s1, t0, t0);
  s2 += (t0 + t1);

  /* O(eps^3) order terms */
  p6 = qd::two_prod(a[0], b[3], q6);
  p7 = qd::two_prod(a[1], b[2], q7);
  p8 = qd::two_prod(a[2], b[1], q8);
  p9 = qd::two_prod(a[3], b[0], q9);

  /* Nine-Two-Sum of q0, s1, q3, q4, q5, p6, p7, p8, p9. */
  q0 = qd::two_sum(q0, q3, q3);
  q4 = qd::two_sum(q4, q5, q5);
  p6 = qd::two_sum(p6, p7, p7);
  p8 = qd::two_sum(p8, p9, p9);
  /* Compute (t0, t1) = (q0, q3) + (q4, q5). */
  t0 = qd::two_sum(q0, q4, t1);
  t1 += (q3 + q5);
  /* Compute (r0, r1) = (p6, p7) + (p8, p9). */
  r0 = qd::two_sum(p6, p8, r1);
  r1 += (p7 + p9);
  /* Compute (q3, q4) = (t0, t1) + (r0, r1). */
  q3 = qd::two_sum(t0, r0, q4);
  q4 += (t1 + r1);
  /* Compute (t0, t1) = (q3, q4) + s1. */
  t0 = qd::two_sum(q3, s1, t1);
  t1 += q4;

  /* O(eps^4) terms -- Nine-One-Sum */
  t1 += a[1] * b[3] + a[2] * b[2] + a[3] * b[1] + q6 + q7 + q8 + q9 + s2;

  qd::renorm(p0, p1, s0, t0, t1);
  return qd_real(p0, p1, s0, t0);
}

// 2025-12-24(Wed) T.Kouya
// Branch free algorithm
// void Mul4(const double x[4], const double y[4], double z[4]) {
//static inline void c_qd_mul_bf(const double *a, const double *b, double *c)
inline qd_real qd_real::bf_mul(const qd_real &a, const qd_real &b)
{
	double a0, b0, c0, d0, e0, f0, g0, h0, i0, j0, k0, l0, m0, n0, o0, p0;
	double a1, b1, c1, d1, e1, f1, g1, h1, i1, j1, k1, l1, m1, n1;
	double a2, b2, c2, d2, e2, f2, g2, h2, i2, j2, k2, l2, m2;
	double a3, b3, c3, d3, e3, f3, g3, h3;
	double a4, b4, c4, d4, e4, f4;
	double a5, b5, c5, d5;
	double a6, b6, c6, d6;
	double a7, b7, c7, d7;
	double a8, b8, c8, d8;
	double a9, b9, c9, d9;
	double a10, b10, c10, d10;

  a0 = qd::two_prod(a[0], b[0], b0);
  c0 = qd::two_prod(a[0], b[1], e0);
  d0 = qd::two_prod(a[1], b[0], f0);
  g0 = qd::two_prod(a[0], b[2], j0);
  h0 = qd::two_prod(a[1], b[1], k0);
  i0 = qd::two_prod(a[2], b[0], l0);
  m0 = a[0] * b[3];
  n0 = a[1] * b[2];
  o0 = a[2] * b[1];
  p0 = a[3] * b[0];
  c1 = qd::two_sum(c0, d0, d1);
  e1 = qd::two_sum(e0, f0, f1);
  g1 = qd::two_sum(g0, i0, i1);
  j1 = j0 + l0;
  m1 = m0 + p0;
  n1 = n0 + o0;
  b2 = qd::two_sum(b0, c1, c2);
  e2 = qd::two_sum(e1, h0, h2);
  f2 = f1 + j1;
  i2 = i1 + k0;
  m2 = m1 + n1;
  a3 = qd::quick_two_sum(a0, b2, b3);
  c3 = qd::quick_two_sum(c2, d1, d3);
  e3 = qd::two_sum(e2, g1, g3);
  f3 = f2 + m2;
  h3 = h2 + i2;
  c4 = qd::two_sum(c3, e3, e4);
  d4 = d3 + h3;
  f4 = f3 + g3;
  d5 = d4 + e4;
  c6 = qd::two_sum(c4, d5, d6);
  b7 = qd::two_sum(b3, c6, c7);
  d7 = d6 + f4;
  a8 = qd::quick_two_sum(a3, b7, b8);
  c8 = qd::two_sum(c7, d7, d8);
  b9 = qd::two_sum(b8, c8, c9);
  c10 = qd::quick_two_sum(c9, d8, d10);

	//c[0] = a8;
	//c[1] = b9;
	//c[2] = c10;
	//c[3] = d10;
  return qd_real(a8, b9, c10, d10);
}

inline qd_real operator*(const qd_real &a, const qd_real &b) {
#ifdef QD_SLOPPY_MUL
  //return qd_real::sloppy_mul(a, b);
  return qd_real::bf_mul(a, b);
#else
  return qd_real::accurate_mul(a, b);
#endif
}

/* QW-FMA  z = a * b + c   (176 flops)
   Machine-proved with FPANVerifier + z3 5.0.0 (ACS2026 formulation):
     error bound      |z-(ab+c)| <= 812 u^4 (|ab|+|c|)
     every FastTwoSum precondition  exp(x) >= exp(y)
     non-overlapping output         z0 |> z1 |> z2 |> z3   (strongly_dominates)
   The precision p is symbolic, so the proofs hold for binary32/64/128 alike.
   Normalization repeats a cascade over adjacent pairs; the pass count is the
   smallest for which the non-overlap is provable (DW 1 / TW 3 / QW 5).  Two
   passes -- the classical choice -- cannot prove it: for QW not even the
   weakest relation p_dominates holds, because TwoSum(P00,c0) may cancel
   completely when a*b ~= -c and the leading word must be rebuilt from below. */
inline qd_real qw_fma(const qd_real &a, const qd_real &b, const qd_real &c) {
  double P00, E00, P01, E01, P10, E10, P02, E02, P11, E11, P20, E20;
  double P03, P12, P21, P30, D, B, r, Et;
  double A1, f1, f2, f3, f4;
  double A2, g1, g2, g3, g4, g5, g6, g7, g8, g9;
  double A3, t1, t2, t3, t4;
  double w0, w1, w2, w3;
  P00 = qd::two_prod(a[0], b[0], E00);
  P01 = qd::two_prod(a[0], b[1], E01);
  P10 = qd::two_prod(a[1], b[0], E10);
  P02 = qd::two_prod(a[0], b[2], E02);
  P11 = qd::two_prod(a[1], b[1], E11);
  P20 = qd::two_prod(a[2], b[0], E20);
  P03 = a[0] * b[3];
  P12 = a[1] * b[2];
  P21 = a[2] * b[1];
  P30 = a[3] * b[0];
  D   = (P03 + P30) + (P12 + P21);
  B   = qd::two_sum(P00, c[0], r);
  A1  = qd::two_sum(P01, P10, f1);
  A1  = qd::two_sum(A1, E00, f2);
  A1  = qd::two_sum(A1, c[1], f3);
  A1  = qd::two_sum(A1, r, f4);
  A2  = qd::two_sum(P02, P20, g1);
  A2  = qd::two_sum(A2, P11, g2);
  Et  = qd::two_sum(E01, E10, g4);
  A2  = qd::two_sum(A2, Et, g3);
  A2  = qd::two_sum(A2, c[2], g5);
  A2  = qd::two_sum(A2, f1, g6);
  A2  = qd::two_sum(A2, f2, g7);
  A2  = qd::two_sum(A2, f3, g8);
  A2  = qd::two_sum(A2, f4, g9);
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
  double z0, z1, z2, z3;
  /* renormalization pass 1/5 */
  w0 = qd::quick_two_sum(B, A1, w1);
  w1 = qd::two_sum(w1, A2, w2);
  w2 = qd::two_sum(w2, A3, w3);
  /* renormalization pass 2/5 */
  w0 = qd::two_sum(w0, w1, w1);
  w1 = qd::two_sum(w1, w2, w2);
  w2 = qd::quick_two_sum(w2, w3, w3);
  /* renormalization pass 3/5 */
  w0 = qd::two_sum(w0, w1, w1);
  w1 = qd::quick_two_sum(w1, w2, w2);
  w2 = qd::quick_two_sum(w2, w3, w3);
  /* renormalization pass 4/5 */
  w0 = qd::quick_two_sum(w0, w1, w1);
  w1 = qd::quick_two_sum(w1, w2, w2);
  w2 = qd::quick_two_sum(w2, w3, w3);
  /* renormalization pass 5/5 */
  z0 = qd::quick_two_sum(w0, w1, w1);
  z1 = qd::quick_two_sum(w1, w2, w2);
  z2 = qd::quick_two_sum(w2, w3, z3);
  return qd_real(z0, z1, z2, z3);
}

/* QW-FMA  z = a * b + c   (176 flops, scalar multiplier)
   Machine-proved with FPANVerifier + z3 5.0.0 (ACS2026 formulation):
     error bound      |z-(ab+c)| <= 812 u^4 (|ab|+|c|)
     every FastTwoSum precondition  exp(x) >= exp(y)
     non-overlapping output         z0 |> z1 |> z2 |> z3   (strongly_dominates)
   The precision p is symbolic, so the proofs hold for binary32/64/128 alike.
   Normalization repeats a cascade over adjacent pairs; the pass count is the
   smallest for which the non-overlap is provable (DW 1 / TW 3 / QW 5).  Two
   passes -- the classical choice -- cannot prove it: for QW not even the
   weakest relation p_dominates holds, because TwoSum(P00,c0) may cancel
   completely when a*b ~= -c and the leading word must be rebuilt from below. */

/* div/sqrt-safe variants: the Newton iterations of division and square root
   receive residuals that are NOT non-overlapping expansions, so the input
   assumptions behind the machine proof do not hold and no FastTwoSum
   precondition can be claimed for any gate.  A FastTwoSum whose precondition
   fails does not even satisfy s+e=a+b -- it is no longer an EFT -- so these
   variants use TwoSum everywhere.  The pass count is the same as the standard
   variant, hence a safe variant is never cheaper than the standard one
   (DW 20 / TW 84 / QW 206 flops). */
inline qd_real qw_fma_safe(const qd_real &a, double b, const qd_real &c) {
  double P00, E00, P01, E01, P10, E10, P02, E02, P11, E11, P20, E20;
  double P03, P12, P21, P30, D, B, r, Et;
  double A1, f1, f2, f3, f4;
  double A2, g1, g2, g3, g4, g5, g6, g7, g8, g9;
  double A3, t1, t2, t3, t4;
  double w0, w1, w2, w3;
  P00 = qd::two_prod(a[0], b, E00);
  P01 = 0.0; E01 = 0.0;
  P10 = qd::two_prod(a[1], b, E10);
  P02 = 0.0; E02 = 0.0;
  P11 = 0.0; E11 = 0.0;
  P20 = qd::two_prod(a[2], b, E20);
  P03 = 0.0;
  P12 = 0.0;
  P21 = 0.0;
  P30 = a[3] * b;
  D   = (P03 + P30) + (P12 + P21);
  B   = qd::two_sum(P00, c[0], r);
  A1  = qd::two_sum(P01, P10, f1);
  A1  = qd::two_sum(A1, E00, f2);
  A1  = qd::two_sum(A1, c[1], f3);
  A1  = qd::two_sum(A1, r, f4);
  A2  = qd::two_sum(P02, P20, g1);
  A2  = qd::two_sum(A2, P11, g2);
  Et  = qd::two_sum(E01, E10, g4);
  A2  = qd::two_sum(A2, Et, g3);
  A2  = qd::two_sum(A2, c[2], g5);
  A2  = qd::two_sum(A2, f1, g6);
  A2  = qd::two_sum(A2, f2, g7);
  A2  = qd::two_sum(A2, f3, g8);
  A2  = qd::two_sum(A2, f4, g9);
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
  double z0, z1, z2, z3;
  /* renormalization pass 1/5 */
  w0 = qd::two_sum(B, A1, w1);
  w1 = qd::two_sum(w1, A2, w2);
  w2 = qd::two_sum(w2, A3, w3);
  /* renormalization pass 2/5 */
  w0 = qd::two_sum(w0, w1, w1);
  w1 = qd::two_sum(w1, w2, w2);
  w2 = qd::two_sum(w2, w3, w3);
  /* renormalization pass 3/5 */
  w0 = qd::two_sum(w0, w1, w1);
  w1 = qd::two_sum(w1, w2, w2);
  w2 = qd::two_sum(w2, w3, w3);
  /* renormalization pass 4/5 */
  w0 = qd::two_sum(w0, w1, w1);
  w1 = qd::two_sum(w1, w2, w2);
  w2 = qd::two_sum(w2, w3, w3);
  /* renormalization pass 5/5 */
  z0 = qd::two_sum(w0, w1, w1);
  z1 = qd::two_sum(w1, w2, w2);
  z2 = qd::two_sum(w2, w3, z3);
  return qd_real(z0, z1, z2, z3);
}

inline qd_real qw_fma(const qd_real &a, double b, const qd_real &c) {
  double P00, E00, P01, E01, P10, E10, P02, E02, P11, E11, P20, E20;
  double P03, P12, P21, P30, D, B, r, Et;
  double A1, f1, f2, f3, f4;
  double A2, g1, g2, g3, g4, g5, g6, g7, g8, g9;
  double A3, t1, t2, t3, t4;
  double w0, w1, w2, w3;
  P00 = qd::two_prod(a[0], b, E00);
  P01 = 0.0; E01 = 0.0;
  P10 = qd::two_prod(a[1], b, E10);
  P02 = 0.0; E02 = 0.0;
  P11 = 0.0; E11 = 0.0;
  P20 = qd::two_prod(a[2], b, E20);
  P03 = 0.0;
  P12 = 0.0;
  P21 = 0.0;
  P30 = a[3] * b;
  D   = (P03 + P30) + (P12 + P21);
  B   = qd::two_sum(P00, c[0], r);
  A1  = qd::two_sum(P01, P10, f1);
  A1  = qd::two_sum(A1, E00, f2);
  A1  = qd::two_sum(A1, c[1], f3);
  A1  = qd::two_sum(A1, r, f4);
  A2  = qd::two_sum(P02, P20, g1);
  A2  = qd::two_sum(A2, P11, g2);
  Et  = qd::two_sum(E01, E10, g4);
  A2  = qd::two_sum(A2, Et, g3);
  A2  = qd::two_sum(A2, c[2], g5);
  A2  = qd::two_sum(A2, f1, g6);
  A2  = qd::two_sum(A2, f2, g7);
  A2  = qd::two_sum(A2, f3, g8);
  A2  = qd::two_sum(A2, f4, g9);
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
  double z0, z1, z2, z3;
  /* renormalization pass 1/5 */
  w0 = qd::quick_two_sum(B, A1, w1);
  w1 = qd::two_sum(w1, A2, w2);
  w2 = qd::two_sum(w2, A3, w3);
  /* renormalization pass 2/5 */
  w0 = qd::two_sum(w0, w1, w1);
  w1 = qd::two_sum(w1, w2, w2);
  w2 = qd::quick_two_sum(w2, w3, w3);
  /* renormalization pass 3/5 */
  w0 = qd::two_sum(w0, w1, w1);
  w1 = qd::quick_two_sum(w1, w2, w2);
  w2 = qd::quick_two_sum(w2, w3, w3);
  /* renormalization pass 4/5 */
  w0 = qd::quick_two_sum(w0, w1, w1);
  w1 = qd::quick_two_sum(w1, w2, w2);
  w2 = qd::quick_two_sum(w2, w3, w3);
  /* renormalization pass 5/5 */
  z0 = qd::quick_two_sum(w0, w1, w1);
  z1 = qd::quick_two_sum(w1, w2, w2);
  z2 = qd::quick_two_sum(w2, w3, z3);
  return qd_real(z0, z1, z2, z3);
}

inline qd_real fma(const qd_real &a, const qd_real &b, const qd_real &c) {
  return qw_fma(a, b, c);
}

inline qd_real fma(const qd_real &a, double b, const qd_real &c) {
  return qw_fma(a, b, c);
}

/* quad-double ^ 2  = (x0 + x1 + x2 + x3) ^ 2
                    = x0 ^ 2 + 2 x0 * x1 + (2 x0 * x2 + x1 ^ 2)
                               + (2 x0 * x3 + 2 x1 * x2)           */
inline qd_real sqr(const qd_real &a) {
  double p0, p1, p2, p3, p4, p5;
  double q0, q1, q2, q3;
  double s0, s1;
  double t0, t1;
  
  p0 = qd::two_sqr(a[0], q0);
  p1 = qd::two_prod(2.0 * a[0], a[1], q1);
  p2 = qd::two_prod(2.0 * a[0], a[2], q2);
  p3 = qd::two_sqr(a[1], q3);

  p1 = qd::two_sum(q0, p1, q0);

  q0 = qd::two_sum(q0, q1, q1);
  p2 = qd::two_sum(p2, p3, p3);

  s0 = qd::two_sum(q0, p2, t0);
  s1 = qd::two_sum(q1, p3, t1);

  s1 = qd::two_sum(s1, t0, t0);
  t0 += t1;

  s1 = qd::quick_two_sum(s1, t0, t0);
  p2 = qd::quick_two_sum(s0, s1, t1);
  p3 = qd::quick_two_sum(t1, t0, q0);

  p4 = 2.0 * a[0] * a[3];
  p5 = 2.0 * a[1] * a[2];

  p4 = qd::two_sum(p4, p5, p5);
  q2 = qd::two_sum(q2, q3, q3);

  t0 = qd::two_sum(p4, q2, t1);
  t1 = t1 + p5 + q3;

  p3 = qd::two_sum(p3, t0, p4);
  p4 = p4 + q0 + t1;

  qd::renorm(p0, p1, p2, p3, p4);
  return qd_real(p0, p1, p2, p3);

}

/********** Self-Multiplication **********/
/* quad-double *= double */
inline qd_real &qd_real::operator*=(double a) {
  *this = (*this * a);
  return *this;
}

/* quad-double *= double-double */
inline qd_real &qd_real::operator*=(const dd_real &a) {
  *this = (*this * a);
  return *this;
}

/* quad-double *= quad-double */
inline qd_real &qd_real::operator*=(const qd_real &a) {
  *this = *this * a;
  return *this;
}

inline qd_real operator/ (const qd_real &a, const dd_real &b) {
#ifdef QD_SLOPPY_DIV
  return qd_real::sloppy_div(a, b);
#else
  return qd_real::accurate_div(a, b);
#endif
}

inline qd_real operator/(const qd_real &a, const qd_real &b) {
#ifdef QD_SLOPPY_DIV
  return qd_real::sloppy_div(a, b);
#elif defined(QD_NO_FMA_DIV)
  return qd_real::accurate_div(a, b);
#else
  return qd_real::fma_div(a, b);
#endif
}

/* double / quad-double */
inline qd_real operator/(double a, const qd_real &b) {
  return qd_real(a) / b;
}

/* double-double / quad-double */
inline qd_real operator/(const dd_real &a, const qd_real &b) {
  return qd_real(a) / b;
}

/********** Self-Divisions **********/
/* quad-double /= double */
inline qd_real &qd_real::operator/=(double a) {
  *this = (*this / a);
  return *this;
}

/* quad-double /= double-double */
inline qd_real &qd_real::operator/=(const dd_real &a) {
  *this = (*this / a);
  return *this;
}

/* quad-double /= quad-double */
inline qd_real &qd_real::operator/=(const qd_real &a) {
  *this = (*this / a);
  return *this;
}


/********** Exponentiation **********/
inline qd_real qd_real::operator^(int n) const {
  return pow(*this, n);
}

/********** Miscellaneous **********/
inline qd_real abs(const qd_real &a) {
  return (a[0] < 0.0) ? -a : a;
}

inline qd_real fabs(const qd_real &a) {
  return abs(a);
}

/* Quick version.  May be off by one when qd is very close
   to the middle of two integers.                         */
inline qd_real quick_nint(const qd_real &a) {
  qd_real r = qd_real(qd::nint(a[0]), qd::nint(a[1]), 
      qd::nint(a[2]), qd::nint(a[3]));
  r.renorm();
  return r;
}

/*********** Assignments ************/
/* quad-double = double */
inline qd_real &qd_real::operator=(double a) {
  x[0] = a;
  x[1] = x[2] = x[3] = 0.0;
  return *this;
}

/* quad-double = double-double */
inline qd_real &qd_real::operator=(const dd_real &a) {
  x[0] = a._hi();
  x[1] = a._lo();
  x[2] = x[3] = 0.0;
  return *this;
}

/********** Equality Comparison **********/
inline bool operator==(const qd_real &a, double b) {
  return (a[0] == b && a[1] == 0.0 && a[2] == 0.0 && a[3] == 0.0);
}

inline bool operator==(double a, const qd_real &b) {
  return (b == a);
}

inline bool operator==(const qd_real &a, const dd_real &b) {
  return (a[0] == b._hi() && a[1] == b._lo() && 
          a[2] == 0.0 && a[3] == 0.0);
}

inline bool operator==(const dd_real &a, const qd_real &b) {
  return (b == a);
}

inline bool operator==(const qd_real &a, const qd_real &b) {
  return (a[0] == b[0] && a[1] == b[1] && 
          a[2] == b[2] && a[3] == b[3]);
}


/********** Less-Than Comparison ***********/
inline bool operator<(const qd_real &a, double b) {
  return (a[0] < b || (a[0] == b && a[1] < 0.0));
}

inline bool operator<(double a, const qd_real &b) {
  return (b > a);
}

inline bool operator<(const qd_real &a, const dd_real &b) {
  return (a[0] < b._hi() || 
          (a[0] == b._hi() && (a[1] < b._lo() ||
                            (a[1] == b._lo() && a[2] < 0.0))));
}

inline bool operator<(const dd_real &a, const qd_real &b) {
  return (b > a);
}

inline bool operator<(const qd_real &a, const qd_real &b) {
  return (a[0] < b[0] ||
          (a[0] == b[0] && (a[1] < b[1] ||
                            (a[1] == b[1] && (a[2] < b[2] ||
                                              (a[2] == b[2] && a[3] < b[3]))))));
}

/********** Greater-Than Comparison ***********/
inline bool operator>(const qd_real &a, double b) {
  return (a[0] > b || (a[0] == b && a[1] > 0.0));
}

inline bool operator>(double a, const qd_real &b) {
  return (b < a);
}

inline bool operator>(const qd_real &a, const dd_real &b) {
  return (a[0] > b._hi() || 
          (a[0] == b._hi() && (a[1] > b._lo() ||
                            (a[1] == b._lo() && a[2] > 0.0))));
}

inline bool operator>(const dd_real &a, const qd_real &b) {
  return (b < a);
}

inline bool operator>(const qd_real &a, const qd_real &b) {
  return (a[0] > b[0] ||
          (a[0] == b[0] && (a[1] > b[1] ||
                            (a[1] == b[1] && (a[2] > b[2] ||
                                              (a[2] == b[2] && a[3] > b[3]))))));
}


/********** Less-Than-Or-Equal-To Comparison **********/
inline bool operator<=(const qd_real &a, double b) {
  return (a[0] < b || (a[0] == b && a[1] <= 0.0));
}

inline bool operator<=(double a, const qd_real &b) {
  return (b >= a);
}

inline bool operator<=(const qd_real &a, const dd_real &b) {
  return (a[0] < b._hi() || 
          (a[0] == b._hi() && (a[1] < b._lo() || 
                            (a[1] == b._lo() && a[2] <= 0.0))));
}

inline bool operator<=(const dd_real &a, const qd_real &b) {
  return (b >= a);
}

inline bool operator<=(const qd_real &a, const qd_real &b) {
  return (a[0] < b[0] || 
          (a[0] == b[0] && (a[1] < b[1] ||
                            (a[1] == b[1] && (a[2] < b[2] ||
                                              (a[2] == b[2] && a[3] <= b[3]))))));
}

/********** Greater-Than-Or-Equal-To Comparison **********/
inline bool operator>=(const qd_real &a, double b) {
  return (a[0] > b || (a[0] == b && a[1] >= 0.0));
}

inline bool operator>=(double a, const qd_real &b) {
  return (b <= a);
}

inline bool operator>=(const qd_real &a, const dd_real &b) {
  return (a[0] > b._hi() || 
          (a[0] == b._hi() && (a[1] > b._lo() || 
                            (a[1] == b._lo() && a[2] >= 0.0))));
}

inline bool operator>=(const dd_real &a, const qd_real &b) {
  return (b <= a);
}

inline bool operator>=(const qd_real &a, const qd_real &b) {
  return (a[0] > b[0] || 
          (a[0] == b[0] && (a[1] > b[1] ||
                            (a[1] == b[1] && (a[2] > b[2] ||
                                              (a[2] == b[2] && a[3] >= b[3]))))));
}



/********** Not-Equal-To Comparison **********/
inline bool operator!=(const qd_real &a, double b) {
  return !(a == b);
}

inline bool operator!=(double a, const qd_real &b) {
  return !(a == b);
}

inline bool operator!=(const qd_real &a, const dd_real &b) {
  return !(a == b);
}

inline bool operator!=(const dd_real &a, const qd_real &b) {
  return !(a == b);
}

inline bool operator!=(const qd_real &a, const qd_real &b) {
  return !(a == b);
}



inline qd_real aint(const qd_real &a) {
  return (a[0] >= 0) ? floor(a) : ceil(a);
}

inline bool qd_real::is_zero() const {
  return (x[0] == 0.0);
}

inline bool qd_real::is_one() const {
  return (x[0] == 1.0 && x[1] == 0.0 && x[2] == 0.0 && x[3] == 0.0);
}

inline bool qd_real::is_positive() const {
  return (x[0] > 0.0);
}

inline bool qd_real::is_negative() const {
  return (x[0] < 0.0);
}

inline dd_real to_dd_real(const qd_real &a) {
  return dd_real(a[0], a[1]);
}

inline double to_double(const qd_real &a) {
  return a[0];
}

inline int to_int(const qd_real &a) {
  return static_cast<int>(a[0]);
}

inline qd_real inv(const qd_real &qd) {
  return 1.0 / qd;
}

inline qd_real max(const qd_real &a, const qd_real &b) {
  return (a > b) ? a : b;
}

inline qd_real max(const qd_real &a, const qd_real &b, 
                   const qd_real &c) {
  return (a > b) ? ((a > c) ? a : c) : ((b > c) ? b : c);
}

inline qd_real min(const qd_real &a, const qd_real &b) {
  return (a < b) ? a : b;
}

inline qd_real min(const qd_real &a, const qd_real &b, 
                   const qd_real &c) {
  return (a < b) ? ((a < c) ? a : c) : ((b < c) ? b : c);
}

/* Random number generator */
inline qd_real qd_real::rand() {
  return qdrand();
}

inline qd_real ldexp(const qd_real &a, int n) {
  return qd_real(std::ldexp(a[0], n), std::ldexp(a[1], n), 
                 std::ldexp(a[2], n), std::ldexp(a[3], n));
}

#endif /* _QD_QD_INLINE_H */
