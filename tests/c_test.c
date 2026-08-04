#include <stdio.h>
#include <math.h>
#include <qd/c_dd.h>
#include <qd/c_td.h>
#include <qd/c_qd.h>

/* Test 1 (dd).  Salamin-Brent quadratically convergent formula for pi. */
static int test_1_dd(void) {

  double a[2], b[2], s[2], p[2], t[2], t2[2];
  double a_new[2], b_new[2], p_old[2], tmp[2];
  double m;
  int r, i;
  const int max_iter = 20;

  puts("Test 1 (dd_real).  (Salamin-Brent quadratically convergent formula for pi)");

  c_dd_copy_d(1.0, a);  /* a = 1.0 */
  c_dd_copy_d(0.5, t);  /* t = 0.5 */
  c_dd_sqrt(t, b);      /* b = sqrt(t) */
  c_dd_copy_d(0.5, s);  /* s = 0.5 */
  m = 1.0;

  c_dd_sqr(a, p);
  c_dd_mul_dd_d(p, 2.0, tmp); c_dd_copy(tmp, p);
  c_dd_div(p, s, tmp);        c_dd_copy(tmp, p);

  printf("  iteration 0: ");
  c_dd_write(p);
  for (i = 1; i <= max_iter; i++) {
    m *= 2.0;

    /* a_new = 0.5 * (a + b) */
    c_dd_add(a, b, a_new);
    c_dd_mul_dd_d(a_new, 0.5, tmp); c_dd_copy(tmp, a_new);

    c_dd_mul(a, b, b_new); /* b_new = a * b */

    /* Compute s = s - m * (a_new^2 - b_new) */
    c_dd_sqr(a_new, t);                /* t = a_new ^ 2 */
    c_dd_sub(t, b_new, tmp); c_dd_copy(tmp, t);  /* t -= b_new */
    c_dd_mul_dd_d(t, m, tmp); c_dd_copy(tmp, t); /* t *= m */
    c_dd_sub(s, t, tmp); c_dd_copy(tmp, s);      /* s -= t */

    c_dd_copy(a_new, a);
    c_dd_sqrt(b_new, b);
    c_dd_copy(p, p_old);

    /* Compute  p = 2.0 * a^2 / s */
    c_dd_sqr(a, p);
    c_dd_mul_dd_d(p, 2.0, tmp); c_dd_copy(tmp, p);
    c_dd_div(p, s, tmp);        c_dd_copy(tmp, p);

    /* Test for convergence by looking at |p - p_old|. */
    c_dd_sub(p, p_old, t);
    c_dd_abs(t, t2);
    c_dd_comp_dd_d(t2, 1e-30, &r);
    if (r < 0) break;

    printf("  iteration %1d: ", i);
    c_dd_write(p);
  }

  c_dd_pi(p);   /* p = pi */
  printf("          _pi: ");
  c_dd_write(p);
  printf("        error: %.5e = %g eps\n", t2[0], t2[0] / ldexp(1.0, -103));

  return 0;
}

/* Test 1 (td).  Salamin-Brent quadratically convergent formula for pi. */
static int test_1_td(void) {

  double a[3], b[3], s[3], p[3], t[3], t2[3];
  double a_new[3], b_new[3], p_old[3];
  double m;
  int r, i;
  const int max_iter = 20;

  puts("Test 1 (td_real).  (Salamin-Brent quadratically convergent formula for pi)");

  c_td_copy_d(1.0, a);  /* a = 1.0 */
  c_td_copy_d(0.5, t);  /* t = 0.5 */
  c_td_sqrt(t, b);      /* b = sqrt(t) */
  c_td_copy_d(0.5, s);  /* s = 0.5 */
  m = 1.0;

  c_td_sqr(a, p);
  c_td_selfmul_d(2.0, p);
  c_td_selfdiv(s, p);

  printf("  iteration 0: ");
  c_td_write(p);
  for (i = 1; i <= max_iter; i++) {
    m *= 2.0;

    /* a_new = 0.5 * (a + b) */
    c_td_add(a, b, a_new);
    c_td_selfmul_d(0.5, a_new);

    c_td_mul(a, b, b_new); /* b_new = a * b */

    /* Compute s = s - m * (a_new^2 - b_new) */
    c_td_sqr(a_new, t);       /* t = a_new ^ 2 */
    c_td_selfsub(b_new, t);   /* t -= b_new */
    c_td_selfmul_d(m, t);     /* t *= m */
    c_td_selfsub(t, s);       /* s -= t */

    c_td_copy(a_new, a);
    c_td_sqrt(b_new, b);
    c_td_copy(p, p_old);

    /* Compute  p = 2.0 * a^2 / s */
    c_td_sqr(a, p);
    c_td_selfmul_d(2.0, p);
    c_td_selfdiv(s, p);

    /* Test for convergence by looking at |p - p_old|. */
    c_td_sub(p, p_old, t);
    c_td_abs(t, t2);
    c_td_comp_td_d(t2, 1e-46, &r);
    if (r < 0) break;

    printf("  iteration %1d: ", i);
    c_td_write(p);
  }

  c_td_pi(p);   /* p = pi */
  printf("          _pi: ");
  c_td_write(p);
  printf("        error: %.5e = %g eps\n", t2[0], t2[0] / ldexp(1.0, -156));

  return 0;
}

/* Test 1 (qd).  Salamin-Brent quadratically convergent formula for pi. */
static int test_1_qd(void) {

  double a[4], b[4], s[4], p[4], t[4], t2[4];
  double a_new[4], b_new[4], p_old[4];
  double m;
  int r, i;
  const int max_iter = 20;

  puts("Test 1 (qd_real).  (Salamin-Brent quadratically convergent formula for pi)");

  c_qd_copy_d(1.0, a);  /* a = 1.0 */
  c_qd_copy_d(0.5, t);  /* t = 0.5 */
  c_qd_sqrt(t, b);      /* b = sqrt(t) */
  c_qd_copy_d(0.5, s);  /* s = 0.5 */
  m = 1.0;

  c_qd_sqr(a, p);
  c_qd_selfmul_d(2.0, p);
  c_qd_selfdiv(s, p);

  printf("  iteration 0: ");
  c_qd_write(p);
  for (i = 1; i <= max_iter; i++) {
    m *= 2.0;

    /* a_new = 0.5 * (a + b) */
    c_qd_add(a, b, a_new);
    c_qd_selfmul_d(0.5, a_new);

    c_qd_mul(a, b, b_new); /* b_new = a * b */

    /* Compute s = s - m * (a_new^2 - b_new) */
    c_qd_sqr(a_new, t);       /* t = a_new ^ 2 */
    c_qd_selfsub(b_new, t);   /* t -= b_new */
    c_qd_selfmul_d(m, t);     /* t *= m */
    c_qd_selfsub(t, s);       /* s -= t */

    c_qd_copy(a_new, a);
    c_qd_sqrt(b_new, b);
    c_qd_copy(p, p_old);

    /* Compute  p = 2.0 * a^2 / s */
    c_qd_sqr(a, p);
    c_qd_selfmul_d(2.0, p);
    c_qd_selfdiv(s, p);

    /* Test for convergence by looking at |p - p_old|. */
    c_qd_sub(p, p_old, t);
    c_qd_abs(t, t2);
    c_qd_comp_qd_d(t2, 1e-60, &r);
    if (r < 0) break;

    printf("  iteration %1d: ", i);
    c_qd_write(p);
  }

  c_qd_pi(p);   /* p = pi */
  printf("          _pi: ");
  c_qd_write(p);
  printf("        error: %.5e = %g eps\n", t2[0], t2[0] / ldexp(1.0, -209));

  return 0;
}

/* Test 2.  The C bindings of the branch-free multi-word fused
   multiply-add must agree with a separate multiply and add to within a
   few units in the last place.  */
static int test_2_fma(void) {
  double a[4], b[4], c[4], f[4], g[4], t[4];
  int i, bad = 0;
  double err, tol;

  puts("Test 2.  (branch-free fused multiply-add: c_dd_fma / c_td_fma / c_qd_fma)");

  for (i = 0; i < 4; i++) { a[i] = 0.0; b[i] = 0.0; c[i] = 0.0; }

  /* dd */
  c_dd_copy_d(3.0, a);  c_dd_sqrt(a, a);       /* a = sqrt(3)  */
  c_dd_copy_d(5.0, b);  c_dd_sqrt(b, b);       /* b = sqrt(5)  */
  c_dd_pi(c);
  c_dd_fma(a, b, c, f);
  c_dd_mul(a, b, t);  c_dd_add(t, c, g);
  c_dd_sub(f, g, t);  c_dd_abs(t, t);
  err = t[0] / f[0];
  tol = 8.0 * ldexp(1.0, -104);
  printf("  c_dd_fma vs mul+add: %.3e  (tol %.3e)  %s\n",
         err, tol, (err <= tol) ? "ok" : "FAILED");
  if (err > tol) bad++;

  /* td */
  c_td_copy_d(3.0, a);  c_td_sqrt(a, a);
  c_td_copy_d(5.0, b);  c_td_sqrt(b, b);
  c_td_pi(c);
  c_td_fma(a, b, c, f);
  c_td_mul(a, b, t);  c_td_add(t, c, g);
  c_td_sub(f, g, t);  c_td_abs(t, t);
  err = t[0] / f[0];
  tol = 8.0 * ldexp(1.0, -156);
  printf("  c_td_fma vs mul+add: %.3e  (tol %.3e)  %s\n",
         err, tol, (err <= tol) ? "ok" : "FAILED");
  if (err > tol) bad++;

  /* qd */
  c_qd_copy_d(3.0, a);  c_qd_sqrt(a, a);
  c_qd_copy_d(5.0, b);  c_qd_sqrt(b, b);
  c_qd_pi(c);
  c_qd_fma(a, b, c, f);
  c_qd_mul(a, b, t);  c_qd_add(t, c, g);
  c_qd_sub(f, g, t);  c_qd_abs(t, t);
  err = t[0] / f[0];
  tol = 8.0 * ldexp(1.0, -209);
  printf("  c_qd_fma vs mul+add: %.3e  (tol %.3e)  %s\n",
         err, tol, (err <= tol) ? "ok" : "FAILED");
  if (err > tol) bad++;

  return bad;
}

int main(void) {
  int bad;
  fpu_fix_start(NULL);
  test_1_dd();
  putchar('\n');
  test_1_td();
  putchar('\n');
  test_1_qd();
  putchar('\n');
  bad = test_2_fma();
  return (bad == 0) ? 0 : 1;
}
