/*
 * tests/fma_bench.cpp
 *
 * Benchmark for the branch-free multi-word fused multiply-add
 * (dw_fma / tw_fma / qw_fma) and for the division and square root built
 * on top of them.
 *
 *   "old"  the 0.0.2 code path: a separate multiply followed by an
 *          addition, accurate_div / sloppy_div, sqrt_legacy.
 *   "new"  the fma-based code path that 0.0.3 uses by default.
 *
 * Every kernel runs over arrays of random operands, so nothing can be
 * hoisted out of the loop, and every result is stored into a global
 * array, so nothing can be optimized away.  The number of passes is
 * calibrated at run time and the fastest of several runs is reported.
 *
 * The multiply-and-add is measured twice.  The "independent" form
 *
 *     z[i] = a[i] * b[i] + c[i]
 *
 * is throughput bound and shows what the smaller operation count of the
 * fused form buys.  The "chained" form
 *
 *     s = a[i] * b[i] + s
 *
 * is latency bound (it is the inner loop of a dot product written without
 * multiple accumulators), so it is limited by the dependency chain rather
 * than by the operation count.
 *
 * Created for the dtq package extension (2026).
 */
#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "config.h"
#include <qd/dd_real.h>
#include <qd/td_real.h>
#include <qd/qd_real.h>
#include <qd/ds_real.h>
#include <qd/ts_real.h>
#include <qd/qs_real.h>
#include <qd/fpu.h>

#include "tictoc.h"

using std::printf;

/* Small enough that all four working arrays of any type stay in cache. */
#define N 4096
static const int REPEAT = 5;       /* timed runs; the fastest wins */

static int n_faster = 0, n_slower = 0;

/* Operand and result buffers.  External linkage on purpose: the stores
   into the z arrays must be observable so that no kernel is elided. */
dd_real dd_a[N], dd_b[N], dd_c[N], dd_z[N];
td_real td_a[N], td_b[N], td_c[N], td_z[N];
qd_real qd_a[N], qd_b[N], qd_c[N], qd_z[N];
ds_real ds_a[N], ds_b[N], ds_c[N], ds_z[N];
ts_real ts_a[N], ts_b[N], ts_c[N], ts_z[N];
qs_real qs_a[N], qs_b[N], qs_c[N], qs_z[N];

/* Checksum for the chained kernels. */
double checksum = 0.0;

static void keep(const dd_real &a) { checksum += a.x[0]; }
static void keep(const td_real &a) { checksum += a[0]; }
static void keep(const qd_real &a) { checksum += a[0]; }
static void keep(const ds_real &a) { checksum += a._hi(); }
static void keep(const ts_real &a) { checksum += a[0]; }
static void keep(const qs_real &a) { checksum += a[0]; }

/* label, old kernel statement, new kernel statement.
   Both statements use the loop variable i and must store their result. */
#define BENCH(label, OLD, NEW)                                            \
  do {                                                                    \
    tictoc tv;                                                            \
    tic(&tv);                                                             \
    for (int i = 0; i < N; i++) { OLD; }                                  \
    double t1 = toc(&tv);                                                 \
    int np = (t1 > 0.0) ? static_cast<int>(0.15 / t1) : 200;              \
    if (np < 1) np = 1;                                                   \
                                                                          \
    double t_old = 1e30, t_new = 1e30;                                    \
    for (int r = 0; r < REPEAT; r++) {                                    \
      tic(&tv);                                                           \
      for (int p = 0; p < np; p++)                                        \
        for (int i = 0; i < N; i++) { OLD; }                              \
      double t = toc(&tv) / static_cast<double>(np) / N * 1e9;            \
      if (t < t_old) t_old = t;                                           \
    }                                                                     \
    for (int r = 0; r < REPEAT; r++) {                                    \
      tic(&tv);                                                           \
      for (int p = 0; p < np; p++)                                        \
        for (int i = 0; i < N; i++) { NEW; }                              \
      double t = toc(&tv) / static_cast<double>(np) / N * 1e9;            \
      if (t < t_new) t_new = t;                                           \
    }                                                                     \
    printf("  %-22s %9.2f %9.2f   x%5.2f%s\n", label, t_old, t_new,       \
           t_old / t_new, (t_new < t_old) ? "" : "   <-- no gain");       \
    std::fflush(stdout);                                                  \
    if (t_new < t_old) n_faster++; else n_slower++;                       \
  } while (0)

static void header(const char *title) {
  printf("\n%s\n", title);
  printf("  %-22s %9s %9s   %7s\n", "kernel", "old ns", "new ns", "speedup");
}

static void fill(void) {
  for (int i = 0; i < N; i++) {
    dd_a[i] = ddrand(); dd_b[i] = ddrand(); dd_c[i] = ddrand();
    td_a[i] = tdrand(); td_b[i] = tdrand(); td_c[i] = tdrand();
    qd_a[i] = qdrand(); qd_b[i] = qdrand(); qd_c[i] = qdrand();
    ds_a[i] = dsrand(); ds_b[i] = dsrand(); ds_c[i] = dsrand();
    ts_a[i] = tsrand(); ts_b[i] = tsrand(); ts_c[i] = tsrand();
    qs_a[i] = qsrand(); qs_b[i] = qsrand(); qs_c[i] = qsrand();
  }
}

int main() {
  unsigned int old_cw;
  fpu_fix_start(&old_cw);
  std::srand(20260804);
  fill();

  printf("dtq -- branch-free multi-word FMA benchmark\n");
  printf("%d operands per pass, best of %d runs, ns per operation\n",
         N, REPEAT);

  /* ------- multiply-and-add, independent (throughput bound) ------- */
  header("multiply-and-add:  z = a*b + c   vs   z = dw/tw/qw_fma(a,b,c)");

  BENCH("dd_real  (dw_fma)", dd_z[i] = dd_a[i] * dd_b[i] + dd_c[i],
                             dd_z[i] = dw_fma(dd_a[i], dd_b[i], dd_c[i]));
  BENCH("td_real  (tw_fma)", td_z[i] = td_a[i] * td_b[i] + td_c[i],
                             td_z[i] = tw_fma(td_a[i], td_b[i], td_c[i]));
  BENCH("qd_real  (qw_fma)", qd_z[i] = qd_a[i] * qd_b[i] + qd_c[i],
                             qd_z[i] = qw_fma(qd_a[i], qd_b[i], qd_c[i]));
  BENCH("ds_real  (dw_fma)", ds_z[i] = ds_a[i] * ds_b[i] + ds_c[i],
                             ds_z[i] = dw_fma(ds_a[i], ds_b[i], ds_c[i]));
  BENCH("ts_real  (tw_fma)", ts_z[i] = ts_a[i] * ts_b[i] + ts_c[i],
                             ts_z[i] = tw_fma(ts_a[i], ts_b[i], ts_c[i]));
  BENCH("qs_real  (qw_fma)", qs_z[i] = qs_a[i] * qs_b[i] + qs_c[i],
                             qs_z[i] = qw_fma(qs_a[i], qs_b[i], qs_c[i]));

  /* ------- multiply-and-add, chained (latency bound) ------- */
  header("chained multiply-and-add:  s = a*b + s   (dot-product inner loop)");

  {
    dd_real s(0.0);
    BENCH("dd_real  (dw_fma)", s = dd_a[i] * dd_b[i] + s,
                               s = dw_fma(dd_a[i], dd_b[i], s));
    keep(s);
  }
  {
    td_real s(0.0);
    BENCH("td_real  (tw_fma)", s = td_a[i] * td_b[i] + s,
                               s = tw_fma(td_a[i], td_b[i], s));
    keep(s);
  }
  {
    qd_real s(0.0);
    BENCH("qd_real  (qw_fma)", s = qd_a[i] * qd_b[i] + s,
                               s = qw_fma(qd_a[i], qd_b[i], s));
    keep(s);
  }
  {
    ds_real s(0.0f);
    BENCH("ds_real  (dw_fma)", s = ds_a[i] * ds_b[i] + s,
                               s = dw_fma(ds_a[i], ds_b[i], s));
    keep(s);
  }
  {
    ts_real s(0.0f);
    BENCH("ts_real  (tw_fma)", s = ts_a[i] * ts_b[i] + s,
                               s = tw_fma(ts_a[i], ts_b[i], s));
    keep(s);
  }
  {
    qs_real s(0.0f);
    BENCH("qs_real  (qw_fma)", s = qs_a[i] * qs_b[i] + s,
                               s = qw_fma(qs_a[i], qs_b[i], s));
    keep(s);
  }

  /* ---------------- division ---------------- */
  header("division:  0.0.2 long division   vs   fma_div");

  for (int i = 0; i < N; i++) {   /* keep divisors well away from zero */
    dd_b[i] = abs(dd_b[i]) + 2.0;  td_b[i] = abs(td_b[i]) + 2.0;
    qd_b[i] = abs(qd_b[i]) + 2.0;  ds_b[i] = abs(ds_b[i]) + 2.0f;
    ts_b[i] = abs(ts_b[i]) + 2.0f; qs_b[i] = abs(qs_b[i]) + 2.0f;
  }

  BENCH("dd_real  a/b", dd_z[i] = dd_real::accurate_div(dd_a[i], dd_b[i]),
                        dd_z[i] = dd_real::fma_div(dd_a[i], dd_b[i]));
  BENCH("td_real  a/b", td_z[i] = td_real::accurate_div(td_a[i], td_b[i]),
                        td_z[i] = td_real::fma_div(td_a[i], td_b[i]));
  BENCH("qd_real  a/b", qd_z[i] = qd_real::accurate_div(qd_a[i], qd_b[i]),
                        qd_z[i] = qd_real::fma_div(qd_a[i], qd_b[i]));
  BENCH("ds_real  a/b", ds_z[i] = ds_real::sloppy_div(ds_a[i], ds_b[i]),
                        ds_z[i] = ds_real::fma_div(ds_a[i], ds_b[i]));
  BENCH("ts_real  a/b", ts_z[i] = ts_real::sloppy_div(ts_a[i], ts_b[i]),
                        ts_z[i] = ts_real::fma_div(ts_a[i], ts_b[i]));
  BENCH("qs_real  a/b", qs_z[i] = qs_real::sloppy_div(qs_a[i], qs_b[i]),
                        qs_z[i] = qs_real::fma_div(qs_a[i], qs_b[i]));

  /* ---------------- square root ---------------- */
  header("square root:  sqrt_legacy (0.0.2)   vs   sqrt (fma-based)");

  for (int i = 0; i < N; i++) {
    dd_a[i] = abs(dd_a[i]) + 0.5;  td_a[i] = abs(td_a[i]) + 0.5;
    qd_a[i] = abs(qd_a[i]) + 0.5;  ds_a[i] = abs(ds_a[i]) + 0.5f;
    ts_a[i] = abs(ts_a[i]) + 0.5f; qs_a[i] = abs(qs_a[i]) + 0.5f;
  }

  BENCH("dd_real  sqrt(a)", dd_z[i] = sqrt_legacy(dd_a[i]),
                            dd_z[i] = sqrt(dd_a[i]));
  BENCH("td_real  sqrt(a)", td_z[i] = sqrt_legacy(td_a[i]),
                            td_z[i] = sqrt(td_a[i]));
  BENCH("qd_real  sqrt(a)", qd_z[i] = sqrt_legacy(qd_a[i]),
                            qd_z[i] = sqrt(qd_a[i]));
  BENCH("ds_real  sqrt(a)", ds_z[i] = sqrt_legacy(ds_a[i]),
                            ds_z[i] = sqrt(ds_a[i]));
  BENCH("ts_real  sqrt(a)", ts_z[i] = sqrt_legacy(ts_a[i]),
                            ts_z[i] = sqrt(ts_a[i]));
  BENCH("qs_real  sqrt(a)", qs_z[i] = sqrt_legacy(qs_a[i]),
                            qs_z[i] = sqrt(qs_a[i]));

  printf("\n%d kernel(s) faster, %d kernel(s) not faster.\n",
         n_faster, n_slower);

  fpu_fix_end(&old_cw);
  return 0;
}
