/*
 * include/td_real.h
 *
 * Triple-double precision (>= 156-bit significand) floating point arithmetic
 * package, modeled on the qd_real / dd_real packages by Yozo Hida (LBNL).
 *
 * A td_real is represented by an unevaluated sum of three IEEE doubles
 * x[0] + x[1] + x[2] in non-overlapping form, giving roughly 47 decimal
 * digits of precision.
 *
 * Created for the dtq package extension (2026).
 */
#ifndef _QD_TD_REAL_H
#define _QD_TD_REAL_H

#include <iostream>
#include <string>
#include <limits>
#include <qd/qd_config.h>
#include <qd/dd_real.h>
#include <qd/qd_real.h>

struct QD_API td_real {
  double x[3];    /* The three components. */

  /* Eliminates any zeros in the middle component. */
  void renorm();
  void renorm(double &e);

  td_real(double x0, double x1, double x2);
  explicit td_real(const double *xx);

  static const td_real _2pi;
  static const td_real _pi;
  static const td_real _3pi4;
  static const td_real _pi2;
  static const td_real _pi4;
  static const td_real _e;
  static const td_real _log2;
  static const td_real _log10;
  static const td_real _nan;
  static const td_real _inf;

  static const double _eps;
  static const double _min_normalized;
  static const td_real _max;
  static const td_real _safe_max;
  static const int _ndigits;

  td_real();
  td_real(const char *s);
  td_real(const dd_real &dd);
  td_real(const qd_real &qd);
  td_real(double d);
  td_real(int i);

  double operator[](int i) const;
  double &operator[](int i);

  static void error(const char *msg);

  bool isnan() const;
  bool isfinite() const { return QD_ISFINITE(x[0]); }
  bool isinf() const { return QD_ISINF(x[0]); }

  static td_real ieee_add(const td_real &a, const td_real &b);
  static td_real sloppy_add(const td_real &a, const td_real &b);
  static td_real bf_add(const td_real &a, const td_real &b);

  td_real &operator+=(double a);
  td_real &operator+=(const dd_real &a);
  td_real &operator+=(const td_real &a);

  td_real &operator-=(double a);
  td_real &operator-=(const dd_real &a);
  td_real &operator-=(const td_real &a);

  static td_real sloppy_mul(const td_real &a, const td_real &b);
  static td_real accurate_mul(const td_real &a, const td_real &b);
  static td_real bf_mul(const td_real &a, const td_real &b);

  td_real &operator*=(double a);
  td_real &operator*=(const dd_real &a);
  td_real &operator*=(const td_real &a);

  static td_real sloppy_div(const td_real &a, const dd_real &b);
  static td_real accurate_div(const td_real &a, const dd_real &b);
  static td_real sloppy_div(const td_real &a, const td_real &b);
  static td_real accurate_div(const td_real &a, const td_real &b);

  td_real &operator/=(double a);
  td_real &operator/=(const dd_real &a);
  td_real &operator/=(const td_real &a);

  td_real operator^(int n) const;

  td_real operator-() const;

  td_real &operator=(double a);
  td_real &operator=(const dd_real &a);
  td_real &operator=(const qd_real &a);
  td_real &operator=(const char *s);

  bool is_zero() const;
  bool is_one() const;
  bool is_positive() const;
  bool is_negative() const;

  static td_real rand(void);

  void to_digits(char *s, int &expn, int precision = _ndigits) const;
  void write(char *s, int len, int precision = _ndigits,
      bool showpos = false, bool uppercase = false) const;
  std::string to_string(int precision = _ndigits, int width = 0,
      std::ios_base::fmtflags fmt = static_cast<std::ios_base::fmtflags>(0),
      bool showpos = false, bool uppercase = false, char fill = ' ') const;
  static int read(const char *s, td_real &a);

  /* Debugging methods */
  void dump(const std::string &name = "", std::ostream &os = std::cerr) const;
  void dump_bits(const std::string &name = "",
                 std::ostream &os = std::cerr) const;

  static td_real debug_rand();
};

namespace std {
  template <>
  class numeric_limits<td_real> : public numeric_limits<double> {
  public:
    inline static double epsilon() { return td_real::_eps; }
    inline static double min() { return td_real::_min_normalized; }
    inline static td_real max() { return td_real::_max; }
    inline static td_real safe_max() { return td_real::_safe_max; }
    static const int digits = 156;
    static const int digits10 = 47;
  };
}

QD_API td_real polyeval(const td_real *c, int n, const td_real &x);
QD_API td_real polyroot(const td_real *c, int n,
    const td_real &x0, int max_iter = 48, double thresh = 0.0);

QD_API td_real tdrand(void);
QD_API td_real sqrt(const td_real &a);

QD_API inline bool isnan(const td_real &a) { return a.isnan(); }
QD_API inline bool isfinite(const td_real &a) { return a.isfinite(); }
QD_API inline bool isinf(const td_real &a) { return a.isinf(); }

/* Computes  td * d  where d is known to be a power of 2. */
QD_API td_real mul_pwr2(const td_real &td, double d);

QD_API td_real operator+(const td_real &a, const td_real &b);
QD_API td_real operator+(const dd_real &a, const td_real &b);
QD_API td_real operator+(const td_real &a, const dd_real &b);
QD_API td_real operator+(const td_real &a, double b);
QD_API td_real operator+(double a, const td_real &b);

QD_API td_real operator-(const td_real &a, const td_real &b);
QD_API td_real operator-(const dd_real &a, const td_real &b);
QD_API td_real operator-(const td_real &a, const dd_real &b);
QD_API td_real operator-(const td_real &a, double b);
QD_API td_real operator-(double a, const td_real &b);

QD_API td_real operator*(const td_real &a, const td_real &b);
QD_API td_real operator*(const dd_real &a, const td_real &b);
QD_API td_real operator*(const td_real &a, const dd_real &b);
QD_API td_real operator*(const td_real &a, double b);
QD_API td_real operator*(double a, const td_real &b);

QD_API td_real operator/(const td_real &a, const td_real &b);
QD_API td_real operator/(const dd_real &a, const td_real &b);
QD_API td_real operator/(const td_real &a, const dd_real &b);
QD_API td_real operator/(const td_real &a, double b);
QD_API td_real operator/(double a, const td_real &b);

QD_API td_real sqr(const td_real &a);
QD_API td_real sqrt(const td_real &a);
QD_API td_real pow(const td_real &a, int n);
QD_API td_real pow(const td_real &a, const td_real &b);
QD_API td_real npwr(const td_real &a, int n);

QD_API td_real nroot(const td_real &a, int n);

QD_API td_real rem(const td_real &a, const td_real &b);
QD_API td_real drem(const td_real &a, const td_real &b);
QD_API td_real divrem(const td_real &a, const td_real &b, td_real &r);

dd_real to_dd_real(const td_real &a);
qd_real to_qd_real(const td_real &a);
td_real to_td_real(const qd_real &a);
double  to_double(const td_real &a);
int     to_int(const td_real &a);

QD_API bool operator==(const td_real &a, const td_real &b);
QD_API bool operator==(const td_real &a, const dd_real &b);
QD_API bool operator==(const dd_real &a, const td_real &b);
QD_API bool operator==(double a, const td_real &b);
QD_API bool operator==(const td_real &a, double b);

QD_API bool operator<(const td_real &a, const td_real &b);
QD_API bool operator<(const td_real &a, const dd_real &b);
QD_API bool operator<(const dd_real &a, const td_real &b);
QD_API bool operator<(double a, const td_real &b);
QD_API bool operator<(const td_real &a, double b);

QD_API bool operator>(const td_real &a, const td_real &b);
QD_API bool operator>(const td_real &a, const dd_real &b);
QD_API bool operator>(const dd_real &a, const td_real &b);
QD_API bool operator>(double a, const td_real &b);
QD_API bool operator>(const td_real &a, double b);

QD_API bool operator<=(const td_real &a, const td_real &b);
QD_API bool operator<=(const td_real &a, const dd_real &b);
QD_API bool operator<=(const dd_real &a, const td_real &b);
QD_API bool operator<=(double a, const td_real &b);
QD_API bool operator<=(const td_real &a, double b);

QD_API bool operator>=(const td_real &a, const td_real &b);
QD_API bool operator>=(const td_real &a, const dd_real &b);
QD_API bool operator>=(const dd_real &a, const td_real &b);
QD_API bool operator>=(double a, const td_real &b);
QD_API bool operator>=(const td_real &a, double b);

QD_API bool operator!=(const td_real &a, const td_real &b);
QD_API bool operator!=(const td_real &a, const dd_real &b);
QD_API bool operator!=(const dd_real &a, const td_real &b);
QD_API bool operator!=(double a, const td_real &b);
QD_API bool operator!=(const td_real &a, double b);

QD_API td_real fabs(const td_real &a);
QD_API td_real abs(const td_real &a);    /* same as fabs */

QD_API td_real ldexp(const td_real &a, int n);

QD_API td_real nint(const td_real &a);
QD_API td_real quick_nint(const td_real &a);
QD_API td_real floor(const td_real &a);
QD_API td_real ceil(const td_real &a);
QD_API td_real aint(const td_real &a);

QD_API td_real sin(const td_real &a);
QD_API td_real cos(const td_real &a);
QD_API td_real tan(const td_real &a);
QD_API void sincos(const td_real &a, td_real &s, td_real &c);

QD_API td_real asin(const td_real &a);
QD_API td_real acos(const td_real &a);
QD_API td_real atan(const td_real &a);
QD_API td_real atan2(const td_real &y, const td_real &x);

QD_API td_real exp(const td_real &a);
QD_API td_real log(const td_real &a);
QD_API td_real log10(const td_real &a);

QD_API td_real sinh(const td_real &a);
QD_API td_real cosh(const td_real &a);
QD_API td_real tanh(const td_real &a);
QD_API void sincosh(const td_real &a, td_real &sin_td, td_real &cos_td);

QD_API td_real asinh(const td_real &a);
QD_API td_real acosh(const td_real &a);
QD_API td_real atanh(const td_real &a);

QD_API td_real tdrand(void);

QD_API td_real max(const td_real &a, const td_real &b);
QD_API td_real max(const td_real &a, const td_real &b, const td_real &c);
QD_API td_real min(const td_real &a, const td_real &b);
QD_API td_real min(const td_real &a, const td_real &b, const td_real &c);

QD_API td_real fmod(const td_real &a, const td_real &b);

QD_API std::ostream &operator<<(std::ostream &s, const td_real &a);
QD_API std::istream &operator>>(std::istream &s, td_real &a);
#ifdef QD_INLINE
#include <qd/td_inline.h>
#endif

#endif /* _QD_TD_REAL_H */
