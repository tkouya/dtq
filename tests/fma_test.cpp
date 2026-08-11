/*
 * tests/fma_test.cpp
 *
 * Tests for the branch-free multi-word fused multiply-add routines
 *
 *     dw_fma  (double-word : dd_real, ds_real)
 *     tw_fma  (triple-word : td_real, ts_real)
 *     qw_fma  (quad-word   : qd_real, qs_real)
 *
 * and for the fma-based division and square root built on top of them.
 *
 * The reference value is not another floating-point type: every quantity
 * that appears in  a * b + c  is expanded into an *exact* non-overlapping
 * expansion of doubles (Shewchuk's grow_expansion), the computed result is
 * subtracted from it exactly, and what remains is the true rounding error.
 * Products of two limbs are exact in a double for the float-based types,
 * and exact in two doubles (two_prod) for the double-based types, so the
 * reference is exact for all six classes.
 *
 * Created for the dtq package extension (2026).
 */
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>

#include "config.h"
#include <qd/dd_real.h>
#include <qd/td_real.h>
#include <qd/qd_real.h>
#include <qd/ds_real.h>
#include <qd/ts_real.h>
#include <qd/qs_real.h>
#include <qd/fpu.h>

using std::printf;

static int n_fail = 0;

/* ------------------------------------------------------------------ */
/* Exact expansion arithmetic                                          */
/* ------------------------------------------------------------------ */

#define EXP_MAX 128

/* Add the double b to the non-overlapping expansion e[0..n-1], which is
   kept in increasing order of magnitude.  Returns the new length. */
static int grow(double *e, int n, double b) {
  double Q = b;
  int m = 0;
  for (int i = 0; i < n; i++) {
    double h;
    Q = qd::two_sum(Q, e[i], h);
    if (h != 0.0) e[m++] = h;
  }
  if (Q != 0.0) e[m++] = Q;
  return m;
}

/* |a * b + c - z| computed exactly, where a, b, c and z are given as
   arrays of limbs (already widened to double).  Returns the magnitude of
   the exact residual. */
static double exact_residual(const double *a, int na,
                             const double *b, int nb,
                             const double *c, int nc,
                             const double *z, int nz) {
  double e[EXP_MAX];
  int n = 0;

  for (int i = 0; i < na; i++) {
    for (int j = 0; j < nb; j++) {
      double err;
      double p = qd::two_prod(a[i], b[j], err);
      n = grow(e, n, err);
      n = grow(e, n, p);
    }
  }
  for (int k = 0; k < nc; k++) n = grow(e, n, c[k]);
  for (int k = 0; k < nz; k++) n = grow(e, n, -z[k]);

  /* Sum smallest first; the expansion is non-overlapping, so this is
     accurate to the last bit of the leading term. */
  double s = 0.0;
  for (int i = 0; i < n; i++) s += e[i];
  return std::fabs(s);
}

/* Widen the limbs of any of the six classes into an array of doubles. */
static int limbs(const dd_real &a, double *w) { w[0]=a.x[0]; w[1]=a.x[1]; return 2; }
static int limbs(const td_real &a, double *w) { w[0]=a[0]; w[1]=a[1]; w[2]=a[2]; return 3; }
static int limbs(const qd_real &a, double *w) { w[0]=a[0]; w[1]=a[1]; w[2]=a[2]; w[3]=a[3]; return 4; }
static int limbs(const ds_real &a, double *w) { w[0]=a._hi(); w[1]=a._lo(); return 2; }
static int limbs(const ts_real &a, double *w) { w[0]=a[0]; w[1]=a[1]; w[2]=a[2]; return 3; }
static int limbs(const qs_real &a, double *w) { w[0]=a[0]; w[1]=a[1]; w[2]=a[2]; w[3]=a[3]; return 4; }

/* Exact | a*b + c - z |. */
template <class T>
static double residual(const T &a, const T &b, const T &c, const T &z) {
  double wa[4], wb[4], wc[4], wz[4];
  int na = limbs(a, wa), nb = limbs(b, wb);
  int nc = limbs(c, wc), nz = limbs(z, wz);
  return exact_residual(wa, na, wb, nb, wc, nc, wz, nz);
}

/* Exact | a*b - z |. */
template <class T>
static double residual(const T &a, const T &b, const T &z) {
  double wa[4], wb[4], wz[4];
  int na = limbs(a, wa), nb = limbs(b, wb), nz = limbs(z, wz);
  return exact_residual(wa, na, wb, nb, (double *)0, 0, wz, nz);
}

static double mag(const dd_real &a) { return std::fabs(a.x[0]); }
static double mag(const td_real &a) { return std::fabs(a[0]); }
static double mag(const qd_real &a) { return std::fabs(a[0]); }
static double mag(const ds_real &a) { return std::fabs(a._hi()); }
static double mag(const ts_real &a) { return std::fabs(a[0]); }
static double mag(const qs_real &a) { return std::fabs(a[0]); }

/* ------------------------------------------------------------------ */
/* Reporting                                                           */
/* ------------------------------------------------------------------ */

static void report(const char *what, double worst_ulp, double tol) {
  bool ok = (worst_ulp <= tol);
  printf("  %-34s  max = %7.3f eps   (tol %5.1f)  %s\n",
         what, worst_ulp, tol, ok ? "ok" : "FAILED");
  if (!ok) n_fail++;
}

/* ------------------------------------------------------------------ */
/* Generic tests, instantiated once per class                          */
/* ------------------------------------------------------------------ */

/* a * b + c :  the fused form must agree with the exact value to within a
   few units in the last place of  |a*b| + |c|  (the usual accuracy
   statement for a multiply-and-add; it is not an assertion about
   catastrophic cancellation, which no algorithm can repair). */
template <class T, class S>
static void test_fma(const char *name, T (*rnd)(void), S eps, double tol) {
  double worst = 0.0, worst_naive = 0.0;
  const int N = 20000;

  for (int i = 0; i < N; i++) {
    T a = rnd(), b = rnd(), c = rnd();
    T z = fma(a, b, c);
    T y = a * b + c;

    double scale = mag(a) * mag(b) + mag(c);
    if (scale == 0.0) continue;

    double u = residual(a, b, c, z) / (scale * static_cast<double>(eps));
    double v = residual(a, b, c, y) / (scale * static_cast<double>(eps));
    if (u > worst) worst = u;
    if (v > worst_naive) worst_naive = v;
  }

  char buf[128];
  std::sprintf(buf, "%s  fma(a,b,c)", name);
  report(buf, worst, tol);
  std::sprintf(buf, "%s  a*b+c (unfused)", name);
  report(buf, worst_naive, tol);
}

/* The fused form must be exact whenever the exact result is
   representable.  Limbs are built from small integers scaled by powers of
   two, so a*b+c fits in the target format with room to spare. */
template <class T, class S>
static void test_fma_exact(const char *name, S /*eps*/) {
  int bad = 0;

  for (int i = 0; i < 4000; i++) {
    S a0 = static_cast<S>((std::rand() % 2048) - 1024);
    S b0 = static_cast<S>((std::rand() % 2048) - 1024);
    S c0 = static_cast<S>((std::rand() % 2048) - 1024);
    T a(a0), b(b0), c(c0);
    T z = fma(a, b, c);
    if (residual(a, b, c, z) != 0.0) bad++;
  }

  printf("  %-34s  %s\n", (std::string(name) + "  exact cases").c_str(),
         bad == 0 ? "ok" : "FAILED");
  if (bad) n_fail++;
}

/* Complete cancellation:  with c = -(a*b) the leading words of the fused
   multiply-add cancel almost entirely.  The machine-proved formulation must
   keep its ordinary error bound  |error| <= few * eps * (|a*b| + |c|)  in
   this case too -- this is what lets polyeval run on the plain fma even
   next to a root.  (A FastTwoSum-gated renormalization without the proof
   loses about half of its words here.) */
template <class T, class S>
static void test_fma_cancel(const char *name, T (*rnd)(void), S eps, double tol) {
  double worst = 0.0;
  const int N = 20000;

  for (int i = 0; i < N; i++) {
    T a = rnd(), b = rnd();
    T c = -(a * b);
    T z = fma(a, b, c);

    double scale = mag(a) * mag(b) + mag(c);
    if (scale == 0.0) continue;
    double u = residual(a, b, c, z) / (scale * static_cast<double>(eps));
    if (u > worst) worst = u;
  }

  char buf[128];
  std::sprintf(buf, "%s  cancellation", name);
  report(buf, worst, tol);
}

/* a / b :  b * (a/b) must reproduce a to within a few eps. */
template <class T, class S>
static void test_div(const char *name, T (*rnd)(void), S eps, double tol) {
  double worst = 0.0;

  for (int i = 0; i < 20000; i++) {
    T a = rnd(), b = rnd();
    if (b.is_zero() || a.is_zero()) continue;
    T q = a / b;
    double u = residual(q, b, -a, T(0.0)) / (mag(a) * static_cast<double>(eps));
    if (u > worst) worst = u;
  }

  char buf[128];
  std::sprintf(buf, "%s  a/b", name);
  report(buf, worst, tol);
}

/* sqrt(a) :  sqrt(a)^2 must reproduce a to within a few eps. */
template <class T, class S>
static void test_sqrt(const char *name, T (*rnd)(void), S eps, double tol) {
  double worst = 0.0;

  for (int i = 0; i < 20000; i++) {
    T a = abs(rnd());
    if (a.is_zero()) continue;
    T s = sqrt(a);
    double u = residual(s, s, -a, T(0.0)) / (mag(a) * static_cast<double>(eps));
    if (u > worst) worst = u;
  }

  char buf[128];
  std::sprintf(buf, "%s  sqrt(a)", name);
  report(buf, worst, tol);
}

/* ------------------------------------------------------------------ */

int main() {
  unsigned int old_cw;
  fpu_fix_start(&old_cw);

  std::srand(20260804);

  printf("dw_fma / tw_fma / qw_fma  --  accuracy tests\n");
  printf("(errors are relative to |a*b| + |c|, in units of the type's eps)\n\n");

  printf("double-based classes\n");
  test_fma<dd_real, double>("dd_real (dw_fma)", ddrand, dd_real::_eps, 8.0);
  test_fma<td_real, double>("td_real (tw_fma)", tdrand, td_real::_eps, 8.0);
  test_fma<qd_real, double>("qd_real (qw_fma)", qdrand, qd_real::_eps, 8.0);

  printf("\nfloat-based classes\n");
  test_fma<ds_real, float>("ds_real (dw_fma)", dsrand, ds_real::_eps, 8.0);
  test_fma<ts_real, float>("ts_real (tw_fma)", tsrand, ts_real::_eps, 8.0);
  test_fma<qs_real, float>("qs_real (qw_fma)", qsrand, qs_real::_eps, 8.0);

  printf("\nexactly representable results\n");
  test_fma_exact<dd_real, double>("dd_real (dw_fma)", dd_real::_eps);
  test_fma_exact<td_real, double>("td_real (tw_fma)", td_real::_eps);
  test_fma_exact<qd_real, double>("qd_real (qw_fma)", qd_real::_eps);
  test_fma_exact<ds_real, float>("ds_real (dw_fma)", ds_real::_eps);
  test_fma_exact<ts_real, float>("ts_real (tw_fma)", ts_real::_eps);
  test_fma_exact<qs_real, float>("qs_real (qw_fma)", qs_real::_eps);

  printf("\ncomplete cancellation  (c = -(a*b), errors relative to |a*b| + |c|)\n");
  test_fma_cancel<dd_real, double>("dd_real (dw_fma)", ddrand, dd_real::_eps, 8.0);
  test_fma_cancel<td_real, double>("td_real (tw_fma)", tdrand, td_real::_eps, 8.0);
  test_fma_cancel<qd_real, double>("qd_real (qw_fma)", qdrand, qd_real::_eps, 8.0);
  test_fma_cancel<ds_real, float>("ds_real (dw_fma)", dsrand, ds_real::_eps, 8.0);
  test_fma_cancel<ts_real, float>("ts_real (tw_fma)", tsrand, ts_real::_eps, 8.0);
  test_fma_cancel<qs_real, float>("qs_real (qw_fma)", qsrand, qs_real::_eps, 8.0);

  printf("\ndivision (built on the fused multiply-add)\n");
  test_div<dd_real, double>("dd_real", ddrand, dd_real::_eps, 8.0);
  test_div<td_real, double>("td_real", tdrand, td_real::_eps, 8.0);
  test_div<qd_real, double>("qd_real", qdrand, qd_real::_eps, 8.0);
  test_div<ds_real, float>("ds_real", dsrand, ds_real::_eps, 8.0);
  test_div<ts_real, float>("ts_real", tsrand, ts_real::_eps, 8.0);
  test_div<qs_real, float>("qs_real", qsrand, qs_real::_eps, 8.0);

  printf("\nsquare root (built on the fused multiply-add)\n");
  test_sqrt<dd_real, double>("dd_real", ddrand, dd_real::_eps, 8.0);
  test_sqrt<td_real, double>("td_real", tdrand, td_real::_eps, 8.0);
  test_sqrt<qd_real, double>("qd_real", qdrand, qd_real::_eps, 8.0);
  test_sqrt<ds_real, float>("ds_real", dsrand, ds_real::_eps, 8.0);
  test_sqrt<ts_real, float>("ts_real", tsrand, ts_real::_eps, 8.0);
  test_sqrt<qs_real, float>("qs_real", qsrand, qs_real::_eps, 8.0);

  /* The mixed overloads (multiplier given as a bare double / float) have
     to agree with the full-width ones to within their own rounding. */
  printf("\nmixed-width overloads  fma(a, scalar, c)\n");
  {
    double worst_dd = 0.0, worst_td = 0.0, worst_qd = 0.0;
    for (int i = 0; i < 20000; i++) {
      double s = to_double(ddrand());

      dd_real a = ddrand(), c = ddrand();
      dd_real z = dw_fma(a, s, c);
      double sc = mag(a) * std::fabs(s) + mag(c);
      if (sc > 0.0)
        worst_dd = std::max(worst_dd,
            residual(a, dd_real(s), c, z) / (sc * dd_real::_eps));

      td_real A = tdrand(), C = tdrand();
      td_real Z = tw_fma(A, s, C);
      sc = mag(A) * std::fabs(s) + mag(C);
      if (sc > 0.0)
        worst_td = std::max(worst_td,
            residual(A, td_real(s), C, Z) / (sc * td_real::_eps));

      qd_real X = qdrand(), Y = qdrand();
      qd_real W = qw_fma(X, s, Y);
      sc = mag(X) * std::fabs(s) + mag(Y);
      if (sc > 0.0)
        worst_qd = std::max(worst_qd,
            residual(X, qd_real(s), Y, W) / (sc * qd_real::_eps));
    }
    report("dd_real  dw_fma(a, double, c)", worst_dd, 8.0);
    report("td_real  tw_fma(a, double, c)", worst_td, 8.0);
    report("qd_real  qw_fma(a, double, c)", worst_qd, 8.0);
  }
  {
    double worst_ds = 0.0, worst_ts = 0.0, worst_qs = 0.0;
    for (int i = 0; i < 20000; i++) {
      float s = static_cast<float>(to_double(dsrand()));

      ds_real a = dsrand(), c = dsrand();
      ds_real z = dw_fma(a, s, c);
      double sc = mag(a) * std::fabs(s) + mag(c);
      if (sc > 0.0)
        worst_ds = std::max(worst_ds,
            residual(a, ds_real(s), c, z) / (sc * ds_real::_eps));

      ts_real A = tsrand(), C = tsrand();
      ts_real Z = tw_fma(A, s, C);
      sc = mag(A) * std::fabs(s) + mag(C);
      if (sc > 0.0)
        worst_ts = std::max(worst_ts,
            residual(A, ts_real(s), C, Z) / (sc * ts_real::_eps));

      qs_real X = qsrand(), Y = qsrand();
      qs_real W = qw_fma(X, s, Y);
      sc = mag(X) * std::fabs(s) + mag(Y);
      if (sc > 0.0)
        worst_qs = std::max(worst_qs,
            residual(X, qs_real(s), Y, W) / (sc * qs_real::_eps));
    }
    report("ds_real  dw_fma(a, float, c)", worst_ds, 8.0);
    report("ts_real  tw_fma(a, float, c)", worst_ts, 8.0);
    report("qs_real  qw_fma(a, float, c)", worst_qs, 8.0);
  }

  printf("\n");
  if (n_fail == 0) {
    printf("All fma tests PASSED.\n");
  } else {
    printf("%d fma test(s) FAILED.\n", n_fail);
  }

  fpu_fix_end(&old_cw);
  return (n_fail == 0) ? 0 : 1;
}
