/*
 * src/c_td.cpp
 *
 * Contains C wrapper function for triple-double precision arithmetic.
 * This can be used from fortran code.
 */
#include <cstring>
#include <iostream>

#include "config.h"
#include <qd/td_real.h>
#include <qd/c_td.h>

#define TO_DOUBLE_PTR(a, ptr) ptr[0] = a.x[0]; ptr[1] = a.x[1]; \
                              ptr[2] = a.x[2];

extern "C" {

/* add */
void c_td_add(const double *a, const double *b, double *c) {
  td_real cc;
  cc = td_real(a) + td_real(b);
  TO_DOUBLE_PTR(cc, c);
}
void c_td_add_td_dd(const double *a, const double *b, double *c) {
  td_real cc;
  cc = td_real(a) + dd_real(b);
  TO_DOUBLE_PTR(cc, c);
}
void c_td_add_dd_td(const double *a, const double *b, double *c) {
  td_real cc;
  cc = dd_real(a) + td_real(b);
  TO_DOUBLE_PTR(cc, c);
}
void c_td_add_td_d(const double *a, double b, double *c) {
  td_real cc;
  cc = td_real(a) + b;
  TO_DOUBLE_PTR(cc, c);
}
void c_td_add_d_td(double a, const double *b, double *c) {
  td_real cc;
  cc = a + td_real(b);
  TO_DOUBLE_PTR(cc, c);
}


/* sub */
void c_td_sub(const double *a, const double *b, double *c) {
  td_real cc;
  cc = td_real(a) - td_real(b);
  TO_DOUBLE_PTR(cc, c);
}
void c_td_sub_td_dd(const double *a, const double *b, double *c) {
  td_real cc;
  cc = td_real(a) - dd_real(b);
  TO_DOUBLE_PTR(cc, c);
}
void c_td_sub_dd_td(const double *a, const double *b, double *c) {
  td_real cc;
  cc = dd_real(a) - td_real(b);
  TO_DOUBLE_PTR(cc, c);
}
void c_td_sub_td_d(const double *a, double b, double *c) {
  td_real cc;
  cc = td_real(a) - b;
  TO_DOUBLE_PTR(cc, c);
}
void c_td_sub_d_td(double a, const double *b, double *c) {
  td_real cc;
  cc = a - td_real(b);
  TO_DOUBLE_PTR(cc, c);
}


/* mul */
void c_td_mul(const double *a, const double *b, double *c) {
  td_real cc;
  cc = td_real(a) * td_real(b);
  TO_DOUBLE_PTR(cc, c);
}
void c_td_mul_td_dd(const double *a, const double *b, double *c) {
  td_real cc;
  cc = td_real(a) * dd_real(b);
  TO_DOUBLE_PTR(cc, c);
}
void c_td_mul_dd_td(const double *a, const double *b, double *c) {
  td_real cc;
  cc = dd_real(a) * td_real(b);
  TO_DOUBLE_PTR(cc, c);
}
void c_td_mul_td_d(const double *a, double b, double *c) {
  td_real cc;
  cc = td_real(a) * b;
  TO_DOUBLE_PTR(cc, c);
}
void c_td_mul_d_td(double a, const double *b, double *c) {
  td_real cc;
  cc = a * td_real(b);
  TO_DOUBLE_PTR(cc, c);
}


/* div */
void c_td_div(const double *a, const double *b, double *c) {
  td_real cc;
  cc = td_real(a) / td_real(b);
  TO_DOUBLE_PTR(cc, c);
}
void c_td_div_td_dd(const double *a, const double *b, double *c) {
  td_real cc;
  cc = td_real(a) / dd_real(b);
  TO_DOUBLE_PTR(cc, c);
}
void c_td_div_dd_td(const double *a, const double *b, double *c) {
  td_real cc;
  cc = dd_real(a) / td_real(b);
  TO_DOUBLE_PTR(cc, c);
}
void c_td_div_td_d(const double *a, double b, double *c) {
  td_real cc;
  cc = td_real(a) / b;
  TO_DOUBLE_PTR(cc, c);
}
void c_td_div_d_td(double a, const double *b, double *c) {
  td_real cc;
  cc = a / td_real(b);
  TO_DOUBLE_PTR(cc, c);
}


/* selfadd */
void c_td_selfadd(const double *a, double *b) {
  td_real bb(b);
  bb += td_real(a);
  TO_DOUBLE_PTR(bb, b);
}
void c_td_selfadd_dd(const double *a, double *b) {
  td_real bb(b);
  bb += dd_real(a);
  TO_DOUBLE_PTR(bb, b);
}
void c_td_selfadd_d(double a, double *b) {
  td_real bb(b);
  bb += a;
  TO_DOUBLE_PTR(bb, b);
}


/* selfsub */
void c_td_selfsub(const double *a, double *b) {
  td_real bb(b);
  bb -= td_real(a);
  TO_DOUBLE_PTR(bb, b);
}
void c_td_selfsub_dd(const double *a, double *b) {
  td_real bb(b);
  bb -= dd_real(a);
  TO_DOUBLE_PTR(bb, b);
}
void c_td_selfsub_d(double a, double *b) {
  td_real bb(b);
  bb -= a;
  TO_DOUBLE_PTR(bb, b);
}


/* selfmul */
void c_td_selfmul(const double *a, double *b) {
  td_real bb(b);
  bb *= td_real(a);
  TO_DOUBLE_PTR(bb, b);
}
void c_td_selfmul_dd(const double *a, double *b) {
  td_real bb(b);
  bb *= dd_real(a);
  TO_DOUBLE_PTR(bb, b);
}
void c_td_selfmul_d(double a, double *b) {
  td_real bb(b);
  bb *= a;
  TO_DOUBLE_PTR(bb, b);
}


/* selfdiv */
void c_td_selfdiv(const double *a, double *b) {
  td_real bb(b);
  bb /= td_real(a);
  TO_DOUBLE_PTR(bb, b);
}
void c_td_selfdiv_dd(const double *a, double *b) {
  td_real bb(b);
  bb /= dd_real(a);
  TO_DOUBLE_PTR(bb, b);
}
void c_td_selfdiv_d(double a, double *b) {
  td_real bb(b);
  bb /= a;
  TO_DOUBLE_PTR(bb, b);
}


/* copy */
void c_td_copy(const double *a, double *b) {
  b[0] = a[0];
  b[1] = a[1];
  b[2] = a[2];
}
void c_td_copy_dd(const double *a, double *b) {
  b[0] = a[0];
  b[1] = a[1];
  b[2] = 0.0;
}
void c_td_copy_d(double a, double *b) {
  b[0] = a;
  b[1] = 0.0;
  b[2] = 0.0;
}


void c_td_sqrt(const double *a, double *b) {
  td_real bb;
  bb = sqrt(td_real(a));
  TO_DOUBLE_PTR(bb, b);
}
void c_td_sqr(const double *a, double *b) {
  td_real bb;
  bb = sqr(td_real(a));
  TO_DOUBLE_PTR(bb, b);
}

void c_td_abs(const double *a, double *b) {
  td_real bb;
  bb = abs(td_real(a));
  TO_DOUBLE_PTR(bb, b);
}

void c_td_npwr(const double *a, int n, double *b) {
  td_real bb;
  bb = npwr(td_real(a), n);
  TO_DOUBLE_PTR(bb, b);
}

void c_td_nroot(const double *a, int n, double *b) {
  td_real bb;
  bb = nroot(td_real(a), n);
  TO_DOUBLE_PTR(bb, b);
}

void c_td_nint(const double *a, double *b) {
  td_real bb;
  bb = nint(td_real(a));
  TO_DOUBLE_PTR(bb, b);
}
void c_td_aint(const double *a, double *b) {
  td_real bb;
  bb = aint(td_real(a));
  TO_DOUBLE_PTR(bb, b);
}
void c_td_floor(const double *a, double *b) {
  td_real bb;
  bb = floor(td_real(a));
  TO_DOUBLE_PTR(bb, b);
}
void c_td_ceil(const double *a, double *b) {
  td_real bb;
  bb = ceil(td_real(a));
  TO_DOUBLE_PTR(bb, b);
}

void c_td_log(const double *a, double *b) {
  td_real bb;
  bb = log(td_real(a));
  TO_DOUBLE_PTR(bb, b);
}
void c_td_log10(const double *a, double *b) {
  td_real bb;
  bb = log10(td_real(a));
  TO_DOUBLE_PTR(bb, b);
}
void c_td_exp(const double *a, double *b) {
  td_real bb;
  bb = exp(td_real(a));
  TO_DOUBLE_PTR(bb, b);
}

void c_td_sin(const double *a, double *b) {
  td_real bb;
  bb = sin(td_real(a));
  TO_DOUBLE_PTR(bb, b);
}
void c_td_cos(const double *a, double *b) {
  td_real bb;
  bb = cos(td_real(a));
  TO_DOUBLE_PTR(bb, b);
}
void c_td_tan(const double *a, double *b) {
  td_real bb;
  bb = tan(td_real(a));
  TO_DOUBLE_PTR(bb, b);
}

void c_td_asin(const double *a, double *b) {
  td_real bb;
  bb = asin(td_real(a));
  TO_DOUBLE_PTR(bb, b);
}
void c_td_acos(const double *a, double *b) {
  td_real bb;
  bb = acos(td_real(a));
  TO_DOUBLE_PTR(bb, b);
}
void c_td_atan(const double *a, double *b) {
  td_real bb;
  bb = atan(td_real(a));
  TO_DOUBLE_PTR(bb, b);
}

void c_td_atan2(const double *a, const double *b, double *c) {
  td_real cc;
  cc = atan2(td_real(a), td_real(b));
  TO_DOUBLE_PTR(cc, c);
}

void c_td_sinh(const double *a, double *b) {
  td_real bb;
  bb = sinh(td_real(a));
  TO_DOUBLE_PTR(bb, b);
}
void c_td_cosh(const double *a, double *b) {
  td_real bb;
  bb = cosh(td_real(a));
  TO_DOUBLE_PTR(bb, b);
}
void c_td_tanh(const double *a, double *b) {
  td_real bb;
  bb = tanh(td_real(a));
  TO_DOUBLE_PTR(bb, b);
}

void c_td_asinh(const double *a, double *b) {
  td_real bb;
  bb = asinh(td_real(a));
  TO_DOUBLE_PTR(bb, b);
}
void c_td_acosh(const double *a, double *b) {
  td_real bb;
  bb = acosh(td_real(a));
  TO_DOUBLE_PTR(bb, b);
}
void c_td_atanh(const double *a, double *b) {
  td_real bb;
  bb = atanh(td_real(a));
  TO_DOUBLE_PTR(bb, b);
}

void c_td_sincos(const double *a, double *s, double *c) {
  td_real ss, cc;
  sincos(td_real(a), ss, cc);
  TO_DOUBLE_PTR(cc, c);
  TO_DOUBLE_PTR(ss, s);
}

void c_td_sincosh(const double *a, double *s, double *c) {
  td_real ss, cc;
  sincosh(td_real(a), ss, cc);
  TO_DOUBLE_PTR(cc, c);
  TO_DOUBLE_PTR(ss, s);
}

void c_td_read(const char *s, double *a) {
  td_real aa(s);
  TO_DOUBLE_PTR(aa, a);
}

void c_td_swrite(const double *a, int precision, char *s, int len) {
  td_real(a).write(s, len, precision);
}

void c_td_write(const double *a) {
  std::cout << td_real(a).to_string(td_real::_ndigits) << std::endl;
}

void c_td_neg(const double *a, double *b) {
  b[0] = -a[0];
  b[1] = -a[1];
  b[2] = -a[2];
}

void c_td_rand(double *a) {
  td_real aa;
  aa = tdrand();
  TO_DOUBLE_PTR(aa, a);
}

void c_td_comp(const double *a, const double *b, int *result) {
  td_real aa(a), bb(b);
  if (aa < bb)
    *result = -1;
  else if (aa > bb)
    *result = 1;
  else
    *result = 0;
}

void c_td_comp_td_d(const double *a, double b, int *result) {
  td_real aa(a);
  if (aa < b)
    *result = -1;
  else if (aa > b)
    *result = 1;
  else
    *result = 0;
}

void c_td_comp_d_td(double a, const double *b, int *result) {
  td_real bb(b);
  if (a < bb)
    *result = -1;
  else if (a > bb)
    *result = 1;
  else
    *result = 0;
}

void c_td_pi(double *a) {
  TO_DOUBLE_PTR(td_real::_pi, a);
}

}
