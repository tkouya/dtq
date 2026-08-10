/*
 * elem_bench.cpp
 *
 * dtq-0.0.3 elementary function benchmark:
 *   [1] Exception / special-value handling checks
 *   [2] Accuracy of sqrt / exp / log / trig / inverse-trig / hyperbolic
 *       functions, measured against MPFR references
 *   [3] Performance (ns per call) for all six precision types,
 *       with a native double baseline
 *
 * Types covered: dd_real, td_real, qd_real (double-based)
 *                ds_real, ts_real, qs_real (float-based)
 *
 * Build:  g++ -O2 -ffp-contract=off -I../include elem_bench.cpp \
 *             ../src/.libs/libqd.a -lmpfr -lgmp -o elem_bench
 * Run  :  ./elem_bench [n_accuracy_samples]   (default 200)
 */

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <chrono>
#include <random>
#include <string>
#include <vector>
#include <map>
#include <sstream>
#include <iostream>

#include <mpfr.h>

#include <qd/dd_real.h>
#include <qd/td_real.h>
#include <qd/qd_real.h>
#include <qd/ds_real.h>
#include <qd/ts_real.h>
#include <qd/qs_real.h>
#include <qd/fpu.h>

/* ------------------------------------------------------------------ */
/* Type traits                                                        */
/* ------------------------------------------------------------------ */

template <class R> struct RT;

template <> struct RT<dd_real> {
  static constexpr const char *name = "dd_real";
  static constexpr int limbs = 2, limb_bits = 53;
  static constexpr bool full_math = true;   /* has asin/atan/sinh/... */
};
template <> struct RT<td_real> {
  static constexpr const char *name = "td_real";
  static constexpr int limbs = 3, limb_bits = 53;
  static constexpr bool full_math = true;
};
template <> struct RT<qd_real> {
  static constexpr const char *name = "qd_real";
  static constexpr int limbs = 4, limb_bits = 53;
  static constexpr bool full_math = true;
};
template <> struct RT<ds_real> {
  static constexpr const char *name = "ds_real";
  static constexpr int limbs = 2, limb_bits = 24;
  static constexpr bool full_math = false;  /* sqrt/exp/log/log10/sin/cos/tan only */
};
template <> struct RT<ts_real> {
  static constexpr const char *name = "ts_real";
  static constexpr int limbs = 3, limb_bits = 24;
  static constexpr bool full_math = false;
};
template <> struct RT<qs_real> {
  static constexpr const char *name = "qs_real";
  static constexpr int limbs = 4, limb_bits = 24;
  static constexpr bool full_math = false;
};

template <class R> constexpr int prec_bits() {
  return RT<R>::limbs * RT<R>::limb_bits;
}

/* Exact conversion: unevaluated sum of limbs -> MPFR */
template <class R> void to_mpfr(mpfr_t r, const R &x) {
  mpfr_set_zero(r, 1);
  for (int i = 0; i < RT<R>::limbs; i++)
    mpfr_add_d(r, r, static_cast<double>(x.x[i]), MPFR_RNDN);
}

/* ------------------------------------------------------------------ */
/* Random full-precision input in [lo, hi] (all limbs populated)      */
/* ------------------------------------------------------------------ */

template <class R>
R random_input(std::mt19937_64 &rng, double lo, double hi, bool logscale) {
  std::uniform_real_distribution<double> u01(0.0, 1.0);
  std::uniform_real_distribution<double> um(-0.5, 0.5);
  double a;
  if (logscale) {
    a = lo * std::pow(hi / lo, u01(rng));   /* log-scale assumes 0 < lo < hi */
  } else {
    a = lo + (hi - lo) * u01(rng);
  }
  R x(a);
  /* populate lower limbs with noise so the input uses full precision */
  for (int i = 1; i < RT<R>::limbs; i++) {
    double noise = a * std::ldexp(um(rng), -RT<R>::limb_bits * i);
    if constexpr (RT<R>::limb_bits == 24)
      x += static_cast<float>(noise);
    else
      x += noise;
  }
  return x;
}

/* ------------------------------------------------------------------ */
/* [1] Exception / special-value handling                             */
/* ------------------------------------------------------------------ */

struct CerrCapture {
  std::stringstream ss;
  std::streambuf *old;
  CerrCapture()  { old = std::cerr.rdbuf(ss.rdbuf()); }
  ~CerrCapture() { std::cerr.rdbuf(old); }
};

enum Expect { E_NAN, E_PINF, E_NINF, E_ZERO, E_ONE, E_INFO };

static const char *expect_str(Expect e) {
  switch (e) {
    case E_NAN:  return "NaN";
    case E_PINF: return "+Inf";
    case E_NINF: return "-Inf";
    case E_ZERO: return "0";
    case E_ONE:  return "1";
    default:     return "(report)";
  }
}

static int g_exc_total = 0, g_exc_fail = 0;

/* expected stderr message: M_NO = must be silent, M_YES = must warn,
   M_ANY = don't care (value check only) */
static const int M_NO = 0, M_YES = 1, M_ANY = 2;

template <class R, class F>
void exc_check(const char *label, F f, Expect e, int expect_msg) {
  std::string msg;
  R r;
  {
    CerrCapture cap;
    r = f();
    msg = cap.ss.str();
  }
  /* classify result */
  std::string got;
  bool pass = true;
  if (r.isnan()) {
    got = "NaN";
    if (e != E_INFO) pass = (e == E_NAN);
  } else if (r.isinf()) {
    got = (static_cast<double>(r.x[0]) > 0.0) ? "+Inf" : "-Inf";
    if (e != E_INFO) pass = (e == (got == "+Inf" ? E_PINF : E_NINF));
  } else if (r.is_zero()) {
    got = "0";
    if (e != E_INFO) pass = (e == E_ZERO);
  } else if (r.is_one()) {
    got = "1";
    if (e != E_INFO) pass = (e == E_ONE);
  } else {
    std::ostringstream os;
    os.precision(8);
    os << r;
    got = os.str();
    if (e != E_INFO) pass = false;
  }
  bool msg_seen = !msg.empty();
  bool msg_ok = (expect_msg == M_ANY) || ((expect_msg == M_YES) == msg_seen);
  /* trim message for display */
  std::string m = msg;
  while (!m.empty() && (m.back() == '\n' || m.back() == '\r')) m.pop_back();
  if (m.size() > 44) m = m.substr(0, 41) + "...";

  g_exc_total++;
  bool ok = pass && (e == E_INFO || msg_ok);
  if (!ok) g_exc_fail++;
  std::printf("  %-22s expect=%-8s got=%-10s msg=%-3s  %s %s\n",
              label, expect_str(e), got.c_str(),
              msg_seen ? "yes" : "no",
              (e == E_INFO) ? "INFO" : (ok ? "PASS" : "FAIL"),
              m.empty() ? "" : ("[" + m + "]").c_str());
}

/* checks shared by all six types */
template <class R>
void exc_tests_common() {
  std::printf("--- %s ---\n", RT<R>::name);
  exc_check<R>("sqrt(-1)",   [] { return sqrt(R(-1.0)); },       E_NAN,  true);
  exc_check<R>("sqrt(0)",    [] { return sqrt(R(0.0)); },        E_ZERO, false);
  exc_check<R>("sqrt(NaN)",  [] { return sqrt(R::_nan); },       E_NAN,  false);
  exc_check<R>("log(-1)",    [] { return log(R(-1.0)); },        E_NAN,  true);
  exc_check<R>("exp(+1e5)",  [] { return exp(R(1.0e5)); },       E_PINF, false);
  exc_check<R>("exp(-1e5)",  [] { return exp(R(-1.0e5)); },      E_ZERO, false);
  exc_check<R>("exp(0)",     [] { return exp(R(0.0)); },         E_ONE,  false);
  exc_check<R>("exp(NaN)",   [] { return exp(R::_nan); },        E_NAN,  false);
  exc_check<R>("sin(NaN)",   [] { return sin(R::_nan); },        E_NAN,  M_ANY);
  exc_check<R>("exp(+Inf)",  [] { return exp(R::_inf); },        E_PINF, false);
  exc_check<R>("log(+Inf)",  [] { return log(R::_inf); },        E_INFO, false);
  exc_check<R>("nroot(-2,2)",[] { return nroot(R(-2.0), 2); },   E_NAN,  true);
  exc_check<R>("nroot(2,-1)",[] { return nroot(R(2.0), -1); },   E_NAN,  true);
  exc_check<R>("npwr(0,0)",  [] { return npwr(R(0.0), 0); },     E_NAN,  true);
  exc_check<R>("sin(1e60)",  [] { return sin(R(1.0e60)); },      E_INFO, false);
}

/* double-based types: log(0) -> NaN (+message), full inverse/hyperbolic set */
template <class R>
void exc_tests_full() {
  exc_check<R>("log(0)",     [] { return log(R(0.0)); },         E_NAN,  true);
  exc_check<R>("asin(2)",    [] { return asin(R(2.0)); },        E_NAN,  true);
  exc_check<R>("acos(-2)",   [] { return acos(R(-2.0)); },       E_NAN,  true);
  exc_check<R>("atan2(0,0)", [] { return atan2(R(0.0), R(0.0)); }, E_NAN, true);
  exc_check<R>("acosh(0.5)", [] { return acosh(R(0.5)); },       E_NAN,  true);
  exc_check<R>("atanh(1)",   [] { return atanh(R(1.0)); },       E_NAN,  true);
  exc_check<R>("atanh(-2)",  [] { return atanh(R(-2.0)); },      E_NAN,  true);
}

/* float-based types: log(0) -> -Inf (no message; see ds_real.cpp) */
template <class R>
void exc_tests_float() {
  exc_check<R>("log(0)",     [] { return log(R(0.0)); },         E_NINF, false);
}

/* ------------------------------------------------------------------ */
/* [2] Accuracy vs MPFR                                               */
/* ------------------------------------------------------------------ */

typedef int (*mpfr_fun1)(mpfr_ptr, mpfr_srcptr, mpfr_rnd_t);

struct AccRow {
  std::string func;
  double max_ulp, avg_ulp, min_bits, avg_bits, worst_x;
  int anomalies;
};

template <class R, class F>
AccRow acc_test(const char *fname, F f, mpfr_fun1 mf,
                double lo, double hi, bool logscale,
                int nsamples, std::mt19937_64 &rng) {
  const int pb = prec_bits<R>();
  const int wp = 2 * pb + 64;             /* MPFR working precision */
  const double bits_cap = pb + 24.0;

  mpfr_t mx, mref, mcomp, mtmp;
  mpfr_inits2(wp, mx, mref, mcomp, mtmp, (mpfr_ptr) 0);

  AccRow row;
  row.func = fname;
  row.max_ulp = 0.0; row.avg_ulp = 0.0;
  row.min_bits = bits_cap; row.avg_bits = 0.0;
  row.worst_x = 0.0; row.anomalies = 0;

  for (int i = 0; i < nsamples; i++) {
    R x = random_input<R>(rng, lo, hi, logscale);
    R y = f(x);

    to_mpfr(mx, x);
    mf(mref, mx, MPFR_RNDN);

    if (y.isnan() || y.isinf()) { row.anomalies++; continue; }

    to_mpfr(mcomp, y);
    mpfr_sub(mtmp, mcomp, mref, MPFR_RNDN);
    double relerr;
    if (mpfr_zero_p(mref)) {
      relerr = std::fabs(mpfr_get_d(mtmp, MPFR_RNDN));   /* absolute fallback */
    } else {
      mpfr_div(mtmp, mtmp, mref, MPFR_RNDN);
      relerr = std::fabs(mpfr_get_d(mtmp, MPFR_RNDN));
    }
    double bits = (relerr > 0.0) ? -std::log2(relerr) : bits_cap;
    if (bits > bits_cap) bits = bits_cap;
    double ulp = relerr * std::ldexp(1.0, pb);

    if (ulp > row.max_ulp) { row.max_ulp = ulp; row.worst_x = static_cast<double>(x.x[0]); }
    if (bits < row.min_bits) row.min_bits = bits;
    row.avg_ulp += ulp;
    row.avg_bits += bits;
  }
  int n_ok = nsamples - row.anomalies;
  if (n_ok > 0) { row.avg_ulp /= n_ok; row.avg_bits /= n_ok; }

  mpfr_clears(mx, mref, mcomp, mtmp, (mpfr_ptr) 0);
  return row;
}

static void acc_print_header(const char *tname, int pb, int nsamples) {
  std::printf("--- %s (theoretical %d bits, N=%d) ---\n", tname, pb, nsamples);
  std::printf("  %-8s %12s %12s %10s %10s %13s %5s\n",
              "func", "max_err[ulp]", "avg_err[ulp]",
              "min_bits", "avg_bits", "worst_x", "anom");
}

static void acc_print_row(const AccRow &r) {
  std::printf("  %-8s %12.3g %12.3g %10.1f %10.1f %13.6g %5d\n",
              r.func.c_str(), r.max_ulp, r.avg_ulp,
              r.min_bits, r.avg_bits, r.worst_x, r.anomalies);
}

template <class R>
void acc_tests(int nsamples) {
  std::mt19937_64 rng(20260807);
  acc_print_header(RT<R>::name, prec_bits<R>(), nsamples);

  auto run = [&](const char *n, R (*f)(const R &), mpfr_fun1 mf,
                 double lo, double hi, bool ls) {
    acc_print_row(acc_test<R>(n, f, mf, lo, hi, ls, nsamples, rng));
  };

  run("sqrt",  [](const R &a) { return sqrt(a); },  mpfr_sqrt,  1e-3, 1e3, true);
  run("exp",   [](const R &a) { return exp(a); },   mpfr_exp,   -20.0, 20.0, false);
  run("log",   [](const R &a) { return log(a); },   mpfr_log,   1e-6, 1e6, true);
  run("log10", [](const R &a) { return log10(a); }, mpfr_log10, 1e-6, 1e6, true);
  run("sin",   [](const R &a) { return sin(a); },   mpfr_sin,   -6.28, 6.28, false);
  run("cos",   [](const R &a) { return cos(a); },   mpfr_cos,   -6.28, 6.28, false);
  run("tan",   [](const R &a) { return tan(a); },   mpfr_tan,   -1.5, 1.5, false);

  if constexpr (RT<R>::full_math) {
    run("asin",  [](const R &a) { return asin(a); },  mpfr_asin,  -0.999, 0.999, false);
    run("acos",  [](const R &a) { return acos(a); },  mpfr_acos,  -0.999, 0.999, false);
    run("atan",  [](const R &a) { return atan(a); },  mpfr_atan,  -10.0, 10.0, false);
    run("sinh",  [](const R &a) { return sinh(a); },  mpfr_sinh,  -5.0, 5.0, false);
    run("cosh",  [](const R &a) { return cosh(a); },  mpfr_cosh,  -5.0, 5.0, false);
    run("tanh",  [](const R &a) { return tanh(a); },  mpfr_tanh,  -5.0, 5.0, false);
    run("asinh", [](const R &a) { return asinh(a); }, mpfr_asinh, -10.0, 10.0, false);
    run("acosh", [](const R &a) { return acosh(a); }, mpfr_acosh, 1.001, 100.0, false);
    run("atanh", [](const R &a) { return atanh(a); }, mpfr_atanh, -0.999, 0.999, false);
  }
  std::printf("\n");
}

/* ------------------------------------------------------------------ */
/* [3] Performance                                                    */
/* ------------------------------------------------------------------ */

static volatile double g_sink = 0.0;

/* run one pass over the input array via `pass()`, repeat until enough
   wall time has accumulated, and return ns/call (best pass)          */
template <class F>
double time_pass(F pass, size_t ncalls) {
  using clk = std::chrono::steady_clock;
  pass();                                   /* warm up */
  double best = 1e300;
  double total = 0.0;
  int reps = 0;
  while (total < 0.08 || reps < 3) {
    auto t0 = clk::now();
    pass();
    auto t1 = clk::now();
    double sec = std::chrono::duration<double>(t1 - t0).count();
    total += sec;
    reps++;
    double ns = sec * 1e9 / static_cast<double>(ncalls);
    if (ns < best) best = ns;
    if (reps > 2000) break;
  }
  return best;
}

/* perf[func][type] = ns/call */
typedef std::map<std::string, std::map<std::string, double> > PerfTable;

template <class R, class F>
void perf_one(PerfTable &tab, const char *fname, const std::vector<R> &xs, F f) {
  double ns = time_pass([&] {
    double s = 0.0;
    for (size_t i = 0; i < xs.size(); i++)
      s += static_cast<double>(f(xs[i]).x[0]);
    g_sink += s;
  }, xs.size());
  tab[fname][RT<R>::name] = ns;
}

template <class R>
void perf_tests(PerfTable &tab, size_t n) {
  std::mt19937_64 rng(97);
  std::vector<R> pos, mid, one;
  pos.reserve(n); mid.reserve(n); one.reserve(n);
  for (size_t i = 0; i < n; i++) {
    pos.push_back(random_input<R>(rng, 1e-3, 1e3, true));      /* sqrt, log */
    mid.push_back(random_input<R>(rng, -3.0, 3.0, false));     /* exp, trig */
    one.push_back(random_input<R>(rng, -0.999, 0.999, false)); /* inverse trig */
  }
  perf_one<R>(tab, "sqrt", pos, [](const R &a) { return sqrt(a); });
  perf_one<R>(tab, "exp",  mid, [](const R &a) { return exp(a); });
  perf_one<R>(tab, "log",  pos, [](const R &a) { return log(a); });
  perf_one<R>(tab, "sin",  mid, [](const R &a) { return sin(a); });
  perf_one<R>(tab, "cos",  mid, [](const R &a) { return cos(a); });
  perf_one<R>(tab, "tan",  mid, [](const R &a) { return tan(a); });
  if constexpr (RT<R>::full_math) {
    perf_one<R>(tab, "asin", one, [](const R &a) { return asin(a); });
    perf_one<R>(tab, "atan", mid, [](const R &a) { return atan(a); });
    perf_one<R>(tab, "sinh", mid, [](const R &a) { return sinh(a); });
    perf_one<R>(tab, "tanh", mid, [](const R &a) { return tanh(a); });
  }
  perf_one<R>(tab, "add",  mid, [](const R &a) { return a + a; });
  perf_one<R>(tab, "mul",  mid, [](const R &a) { return a * a; });
  perf_one<R>(tab, "div",  mid, [](const R &a) { return R(1.0) / a; });
}

/* native double baseline */
static void perf_double_baseline(PerfTable &tab, size_t n) {
  std::mt19937_64 rng(97);
  std::uniform_real_distribution<double> upos(1e-3, 1e3), umid(-3.0, 3.0),
                                         uone(-0.999, 0.999);
  std::vector<double> pos(n), mid(n), one(n);
  for (size_t i = 0; i < n; i++) { pos[i] = upos(rng); mid[i] = umid(rng); one[i] = uone(rng); }
  auto run = [&](const char *fn, const std::vector<double> &xs, double (*f)(double)) {
    double ns = time_pass([&] {
      double s = 0.0;
      for (size_t i = 0; i < xs.size(); i++) s += f(xs[i]);
      g_sink += s;
    }, xs.size());
    tab[fn]["double"] = ns;
  };
  run("sqrt", pos, [](double a) { return std::sqrt(a); });
  run("exp",  mid, [](double a) { return std::exp(a); });
  run("log",  pos, [](double a) { return std::log(a); });
  run("sin",  mid, [](double a) { return std::sin(a); });
  run("cos",  mid, [](double a) { return std::cos(a); });
  run("tan",  mid, [](double a) { return std::tan(a); });
  run("asin", one, [](double a) { return std::asin(a); });
  run("atan", mid, [](double a) { return std::atan(a); });
  run("sinh", mid, [](double a) { return std::sinh(a); });
  run("tanh", mid, [](double a) { return std::tanh(a); });
  run("add",  mid, [](double a) { return a + a; });
  run("mul",  mid, [](double a) { return a * a; });
  run("div",  mid, [](double a) { return 1.0 / a; });
}

static void perf_print(const PerfTable &tab) {
  static const char *cols[] = { "double", "dd_real", "td_real", "qd_real",
                                "ds_real", "ts_real", "qs_real" };
  static const char *rows[] = { "add", "mul", "div", "sqrt", "exp", "log",
                                "sin", "cos", "tan", "asin", "atan",
                                "sinh", "tanh" };
  std::printf("  %-6s", "func");
  for (const char *c : cols) std::printf(" %9s", c);
  std::printf("   [ns/call]\n");
  for (const char *r : rows) {
    auto it = tab.find(r);
    if (it == tab.end()) continue;
    std::printf("  %-6s", r);
    for (const char *c : cols) {
      auto jt = it->second.find(c);
      if (jt == it->second.end()) std::printf(" %9s", "-");
      else                        std::printf(" %9.1f", jt->second);
    }
    std::printf("\n");
  }
}

/* ------------------------------------------------------------------ */

int main(int argc, char **argv) {
  unsigned int cw;
  fpu_fix_start(&cw);

  int nacc = (argc > 1) ? std::atoi(argv[1]) : 200;
  if (nacc < 10) nacc = 10;
  const size_t nperf = 4096;

  std::printf("==============================================================\n");
  std::printf(" dtq-0.0.3 elementary function benchmark\n");
  std::printf("==============================================================\n\n");

  std::printf("[1] Exception / special-value handling\n");
  std::printf("    (msg = error message emitted on stderr)\n\n");
  exc_tests_common<dd_real>(); exc_tests_full<dd_real>(); std::printf("\n");
  exc_tests_common<td_real>(); exc_tests_full<td_real>(); std::printf("\n");
  exc_tests_common<qd_real>(); exc_tests_full<qd_real>(); std::printf("\n");
  exc_tests_common<ds_real>(); exc_tests_float<ds_real>(); std::printf("\n");
  exc_tests_common<ts_real>(); exc_tests_float<ts_real>(); std::printf("\n");
  exc_tests_common<qs_real>(); exc_tests_float<qs_real>(); std::printf("\n");
  std::printf("  => %d checks, %d failed\n\n", g_exc_total, g_exc_fail);

  std::printf("[2] Accuracy vs MPFR reference\n");
  std::printf("    err[ulp]: relative error scaled by 2^prec_bits\n");
  std::printf("    min/avg_bits: correct bits (higher is better)\n\n");
  acc_tests<dd_real>(nacc);
  acc_tests<td_real>(nacc);
  acc_tests<qd_real>(nacc);
  acc_tests<ds_real>(nacc);
  acc_tests<ts_real>(nacc);
  acc_tests<qs_real>(nacc);

  std::printf("[3] Performance (best-pass ns/call, array size %zu)\n\n", nperf);
  PerfTable tab;
  perf_double_baseline(tab, nperf);
  perf_tests<dd_real>(tab, nperf);
  perf_tests<td_real>(tab, nperf);
  perf_tests<qd_real>(tab, nperf);
  perf_tests<ds_real>(tab, nperf);
  perf_tests<ts_real>(tab, nperf);
  perf_tests<qs_real>(tab, nperf);
  perf_print(tab);

  fpu_fix_end(&cw);
  return (g_exc_fail == 0) ? 0 : 1;
}
