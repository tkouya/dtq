/*
 * include/c_td.h
 *
 * Contains C wrapper function prototypes for triple-double precision
 * arithmetic.
 */
#ifndef _QD_C_TD_H
#define _QD_C_TD_H

#include <qd/qd_config.h>
#include <qd/fpu.h>

#ifdef __cplusplus
extern "C" {
#endif

/* add */
void c_td_add(const double *a, const double *b, double *c);
void c_td_add_td_dd(const double *a, const double *b, double *c);
void c_td_add_dd_td(const double *a, const double *b, double *c);
void c_td_add_td_d(const double *a, double b, double *c);
void c_td_add_d_td(double a, const double *b, double *c);

/* sub */
void c_td_sub(const double *a, const double *b, double *c);
void c_td_sub_td_dd(const double *a, const double *b, double *c);
void c_td_sub_dd_td(const double *a, const double *b, double *c);
void c_td_sub_td_d(const double *a, double b, double *c);
void c_td_sub_d_td(double a, const double *b, double *c);

/* mul */
void c_td_mul(const double *a, const double *b, double *c);

/* branch-free triple-word fused multiply-add:  d = a * b + c */
void c_td_fma(const double *a, const double *b, const double *c, double *d);
void c_td_fma_td_d(const double *a, double b, const double *c, double *d);
void c_td_mul_td_dd(const double *a, const double *b, double *c);
void c_td_mul_dd_td(const double *a, const double *b, double *c);
void c_td_mul_td_d(const double *a, double b, double *c);
void c_td_mul_d_td(double a, const double *b, double *c);

/* div */
void c_td_div(const double *a, const double *b, double *c);
void c_td_div_td_dd(const double *a, const double *b, double *c);
void c_td_div_dd_td(const double *a, const double *b, double *c);
void c_td_div_td_d(const double *a, double b, double *c);
void c_td_div_d_td(double a, const double *b, double *c);

/* selfadd */
void c_td_selfadd(const double *a, double *b);
void c_td_selfadd_dd(const double *a, double *b);
void c_td_selfadd_d(double a, double *b);

/* selfsub */
void c_td_selfsub(const double *a, double *b);
void c_td_selfsub_dd(const double *a, double *b);
void c_td_selfsub_d(double a, double *b);

/* selfmul */
void c_td_selfmul(const double *a, double *b);
void c_td_selfmul_dd(const double *a, double *b);
void c_td_selfmul_d(double a, double *b);

/* selfdiv */
void c_td_selfdiv(const double *a, double *b);
void c_td_selfdiv_dd(const double *a, double *b);
void c_td_selfdiv_d(double a, double *b);

/* copy */
void c_td_copy(const double *a, double *b);
void c_td_copy_dd(const double *a, double *b);
void c_td_copy_d(double a, double *b);

void c_td_sqrt(const double *a, double *b);
void c_td_sqr(const double *a, double *b);

void c_td_abs(const double *a, double *b);

void c_td_npwr(const double *a, int b, double *c);
void c_td_nroot(const double *a, int b, double *c);

void c_td_nint(const double *a, double *b);
void c_td_aint(const double *a, double *b);
void c_td_floor(const double *a, double *b);
void c_td_ceil(const double *a, double *b);

void c_td_exp(const double *a, double *b);
void c_td_log(const double *a, double *b);
void c_td_log10(const double *a, double *b);

void c_td_sin(const double *a, double *b);
void c_td_cos(const double *a, double *b);
void c_td_tan(const double *a, double *b);

void c_td_asin(const double *a, double *b);
void c_td_acos(const double *a, double *b);
void c_td_atan(const double *a, double *b);
void c_td_atan2(const double *a, const double *b, double *c);

void c_td_sinh(const double *a, double *b);
void c_td_cosh(const double *a, double *b);
void c_td_tanh(const double *a, double *b);

void c_td_asinh(const double *a, double *b);
void c_td_acosh(const double *a, double *b);
void c_td_atanh(const double *a, double *b);

void c_td_sincos(const double *a, double *s, double *c);
void c_td_sincosh(const double *a, double *s, double *c);

void c_td_read(const char *s, double *a);
void c_td_swrite(const double *a, int precision, char *s, int len);
void c_td_write(const double *a);
void c_td_neg(const double *a, double *b);
void c_td_rand(double *a);
void c_td_comp(const double *a, const double *b, int *result);
void c_td_comp_td_d(const double *a, double b, int *result);
void c_td_comp_d_td(double a, const double *b, int *result);
void c_td_pi(double *a);

#ifdef __cplusplus
}
#endif

#endif  /* _QD_C_TD_H */
