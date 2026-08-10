/*
 * src/td_real.cpp
 *
 * Contains implementation of non-inlined functions of triple-double
 * package.  Inlined functions are found in td_inline.h.
 *
 * Most transcendental functions are computed by promoting the argument
 * to qd_real precision, calling the qd_real implementation, and
 * truncating the result back to td_real.  This is correct because
 * qd_real has substantially more precision than td_real and so the
 * round-trip introduces error well below the td_real working precision.
 */
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <string>

#include "config.h"
#include <qd/td_real.h>
#include "util.h"

#include <qd/bits.h>

#ifndef QD_INLINE
#include <qd/td_inline.h>
#endif

using std::cout;
using std::cerr;
using std::endl;
using std::istream;
using std::ostream;
using std::ios_base;
using std::string;
using std::setw;

using namespace qd;

void td_real::error(const char *msg) {
  if (msg) { cerr << "ERROR " << msg << endl; }
}

/********** Rounding **********/

td_real nint(const td_real &a) {
  double x0, x1, x2;

  x0 = nint(a[0]);
  x1 = x2 = 0.0;

  if (x0 == a[0]) {
    /* First double is already an integer. */
    x1 = nint(a[1]);

    if (x1 == a[1]) {
      /* Second double is already an integer. */
      x2 = nint(a[2]);
    } else {
      if (std::abs(x1 - a[1]) == 0.5 && a[2] < 0.0) {
        x1 -= 1.0;
      }
    }
  } else {
    /* First double is not an integer. */
    if (std::abs(x0 - a[0]) == 0.5 && a[1] < 0.0) {
      x0 -= 1.0;
    }
  }

  renorm3(x0, x1, x2);
  return td_real(x0, x1, x2);
}

td_real floor(const td_real &a) {
  double x0, x1, x2;
  x1 = x2 = 0.0;
  x0 = std::floor(a[0]);

  if (x0 == a[0]) {
    x1 = std::floor(a[1]);

    if (x1 == a[1]) {
      x2 = std::floor(a[2]);
    }

    renorm3(x0, x1, x2);
    return td_real(x0, x1, x2);
  }

  return td_real(x0, x1, x2);
}

td_real ceil(const td_real &a) {
  double x0, x1, x2;
  x1 = x2 = 0.0;
  x0 = std::ceil(a[0]);

  if (x0 == a[0]) {
    x1 = std::ceil(a[1]);

    if (x1 == a[1]) {
      x2 = std::ceil(a[2]);
    }

    renorm3(x0, x1, x2);
    return td_real(x0, x1, x2);
  }

  return td_real(x0, x1, x2);
}


/********** Divisions **********/
/* triple-double / double */
td_real operator/(const td_real &a, double b) {
  /* Strategy:  long division using doubles for the quotient digits. */
  double t0, t1;
  double q0, q1, q2;
  td_real r;

  q0 = a[0] / b;

  /* Compute the remainder a - q0 * b */
  t0 = two_prod(q0, b, t1);
  r = a - dd_real(t0, t1);

  /* First correction */
  q1 = r[0] / b;
  t0 = two_prod(q1, b, t1);
  r -= dd_real(t0, t1);

  /* Second correction */
  q2 = r[0] / b;

  renorm3(q0, q1, q2);
  return td_real(q0, q1, q2);
}

td_real::td_real(const char *s) {
  if (td_real::read(s, *this)) {
    td_real::error("(td_real::td_real): INPUT ERROR.");
    *this = td_real::_nan;
  }
}

td_real &td_real::operator=(const char *s) {
  if (td_real::read(s, *this)) {
    td_real::error("(td_real::operator=): INPUT ERROR.");
    *this = td_real::_nan;
  }
  return *this;
}

istream &operator>>(istream &s, td_real &td) {
  char str[255];
  s >> str;
  td = td_real(str);
  return s;
}

ostream &operator<<(ostream &os, const td_real &td) {
  bool showpos = (os.flags() & ios_base::showpos) != 0;
  bool uppercase = (os.flags() & ios_base::uppercase) != 0;
  return os << td.to_string(os.precision(), os.width(), os.flags(),
      showpos, uppercase, os.fill());
}

/* Read a triple-double from s. */
int td_real::read(const char *s, td_real &td) {
  const char *p = s;
  char ch;
  int sign = 0;
  int point = -1;
  int nd = 0;
  int e = 0;
  bool done = false;
  td_real r = 0.0;

  while (*p == ' ') p++;

  while (!done && (ch = *p) != '\0') {
    if (ch >= '0' && ch <= '9') {
      int d = ch - '0';
      r *= 10.0;
      r += static_cast<double>(d);
      nd++;
    } else {
      switch (ch) {
      case '.':
        if (point >= 0) return -1;
        point = nd;
        break;
      case '-':
      case '+':
        if (sign != 0 || nd > 0) return -1;
        sign = (ch == '-') ? -1 : 1;
        break;
      case 'E':
      case 'e':
        int nread;
        nread = std::sscanf(p+1, "%d", &e);
        done = true;
        if (nread != 1) return -1;
        break;
      case ' ':
        done = true;
        break;
      default:
        return -1;
      }
    }
    p++;
  }

  if (point >= 0) {
    e -= (nd - point);
  }

  if (e != 0) {
    r *= (td_real(10.0) ^ e);
  }

  td = (sign < 0) ? -r : r;
  return 0;
}

void td_real::to_digits(char *s, int &expn, int precision) const {
  int D = precision + 1;

  td_real r = abs(*this);
  int e;
  int i, d;

  if (x[0] == 0.0) {
    expn = 0;
    for (i = 0; i < precision; i++) s[i] = '0';
    return;
  }

  /* First determine the (approximate) exponent. */
  e = static_cast<int>(std::floor(std::log10(std::abs(x[0]))));

  if (e < -300) {
    r *= td_real(10.0) ^ 300;
    r /= td_real(10.0) ^ (e + 300);
  } else if (e > 300) {
    r = ldexp(r, -53);
    r /= td_real(10.0) ^ e;
    r = ldexp(r, 53);
  } else {
    r /= td_real(10.0) ^ e;
  }

  if (r >= 10.0) {
    r /= 10.0;
    e++;
  } else if (r < 1.0) {
    r *= 10.0;
    e--;
  }

  if (r >= 10.0 || r < 1.0) {
    td_real::error("(td_real::to_digits): can't compute exponent.");
    return;
  }

  /* Extract the digits */
  for (i = 0; i < D; i++) {
    d = static_cast<int>(r[0]);
    r -= d;
    r *= 10.0;

    s[i] = static_cast<char>(d + '0');
  }

  /* Fix out of range digits. */
  for (i = D-1; i > 0; i--) {
    if (s[i] < '0') {
      s[i-1]--;
      s[i] += 10;
    } else if (s[i] > '9') {
      s[i-1]++;
      s[i] -= 10;
    }
  }

  if (s[0] <= '0') {
    td_real::error("(td_real::to_digits): non-positive leading digit.");
    return;
  }

  /* Round, handle carry */
  if (s[D-1] >= '5') {
    s[D-2]++;

    i = D-2;
    while (i > 0 && s[i] > '9') {
      s[i] -= 10;
      s[--i]++;
    }
  }

  if (s[0] > '9') {
    e++;
    for (i = precision; i >= 2; i--) s[i] = s[i-1];
    s[0] = '1';
    s[1] = '0';
  }

  s[precision] = 0;
  expn = e;
}

void td_real::write(char *s, int len, int precision,
    bool showpos, bool uppercase) const {
  string str = to_string(precision, 0, ios_base::scientific, showpos, uppercase);
  strncpy(s, str.c_str(), len-1);
  s[len-1] = 0;
}

static void round_string_td(char *s, int precision, int *offset) {
  int i;
  int D = precision;

  if (s[D-1] >= '5') {
    s[D-2]++;

    i = D-2;
    while (i > 0 && s[i] > '9') {
      s[i] -= 10;
      s[--i]++;
    }
  }

  if (s[0] > '9') {
    for (i = precision; i >= 2; i--) s[i] = s[i-1];
    s[0] = '1';
    s[1] = '0';

    (*offset)++;
    precision++;
  }

  s[precision] = 0;
}


string td_real::to_string(int precision, int width, ios_base::fmtflags fmt,
    bool showpos, bool uppercase, char fill) const {
  string s;
  bool fixed = (fmt & ios_base::fixed) != 0;
  bool sgn = true;
  int i, e = 0;

  if (isinf()) {
    if (*this < 0.0)
      s += '-';
    else if (showpos)
      s += '+';
    else
      sgn = false;
    s += uppercase ? "INF" : "inf";
  } else if (isnan()) {
    s = uppercase ? "NAN" : "nan";
    sgn = false;
  } else {
    if (*this < 0.0)
      s += '-';
    else if (showpos)
      s += '+';
    else
      sgn = false;

    if (*this == 0.0) {
      s += '0';
      if (precision > 0) {
        s += '.';
        s.append(precision, '0');
      }
    } else {
      int off = (fixed ? (1 + to_int(floor(log10(abs(*this))))) : 1);
      int d = precision + off;

      int d_with_extra = d;
      if (fixed)
        d_with_extra = std::max(90, d);

      if (fixed && (precision == 0) && (abs(*this) < 1.0)) {
        if (abs(*this) >= 0.5)
          s += '1';
        else
          s += '0';
        return s;
      }

      if (fixed && d <= 0) {
        s += '0';
        if (precision > 0) {
          s += '.';
          s.append(precision, '0');
        }
      } else {
        char *t;
        int j;

        if (fixed) {
          t = new char[d_with_extra+1];
          to_digits(t, e, d_with_extra);
        } else {
          t = new char[d+1];
          to_digits(t, e, d);
        }

        if (fixed) {
          round_string_td(t, d + 1, &off);

          if (off > 0) {
            for (i = 0; i < off; i++) s += t[i];
            if (precision > 0) {
              s += '.';
              for (j = 0; j < precision; j++, i++) s += t[i];
            }
          } else {
            s += "0.";
            if (off < 0) s.append(-off, '0');
            for (i = 0; i < d; i++) s += t[i];
          }
        } else {
          s += t[0];
          if (precision > 0) s += '.';

          for (i = 1; i <= precision; i++)
            s += t[i];
        }
        delete [] t;
      }
    }

    if (fixed && (precision > 0)) {
      double from_string = atof(s.c_str());

      if (fabs(from_string / this->x[0]) > 3.0) {

        for (i = 1; i < (int)s.length(); i++) {
          if (s[i] == '.') {
            s[i] = s[i-1];
            s[i-1] = '.';
            break;
          }
        }

        from_string = atof(s.c_str());
        if (fabs(from_string / this->x[0]) > 3.0) {
          td_real::error("Re-rounding unsuccessful in large number fixed point trap.");
        }
      }
    }

    if (!fixed) {
      s += uppercase ? 'E' : 'e';
      append_expn(s, e);
    }
  }

  /* Fill in the blanks */
  int len = s.length();
  if (len < width) {
    int delta = width - len;
    if (fmt & ios_base::internal) {
      if (sgn)
        s.insert(static_cast<string::size_type>(1), delta, fill);
      else
        s.insert(static_cast<string::size_type>(0), delta, fill);
    } else if (fmt & ios_base::left) {
      s.append(delta, fill);
    } else {
      s.insert(static_cast<string::size_type>(0), delta, fill);
    }
  }

  return s;
}

/* Computes  td^n, where n is an integer. */
td_real pow(const td_real &a, int n) {
  if (n == 0)
    return 1.0;

  td_real r = a;
  td_real s = 1.0;
  int N = std::abs(n);

  if (N > 1) {
    /* Binary exponentiation. */
    while (N > 0) {
      if (N % 2 == 1) {
        s *= r;
      }
      N /= 2;
      if (N > 0)
        r = sqr(r);
    }
  } else {
    s = r;
  }

  if (n < 0)
    return (1.0 / s);

  return s;
}

td_real pow(const td_real &a, const td_real &b) {
  /* Edge cases follow MPFR/IEEE pow semantics so that callers probing an
     endpoint singularity (t^p with t -> 0) get finite answers instead of
     nan from exp(b * log(0)).                                            */
  if (a.is_zero()) {
    if (b.is_zero())
      return 1.0;
    return (b.x[0] > 0.0) ? td_real(0.0) : td_real::_inf;
  }
  return exp(b * log(a));
}

td_real npwr(const td_real &a, int n) {
  return pow(a, n);
}

/* Debugging routines */
void td_real::dump_bits(const string &name, std::ostream &os) const {
  string::size_type len = name.length();
  if (len > 0) {
    os << name << " = ";
    len += 3;
  }
  os << "[ ";
  len += 2;
  for (int j = 0; j < 3; j++) {
    if (j > 0) for (string::size_type i = 0; i < len; i++) os << ' ';
    print_double_info(os, x[j]);
    if (j < 2)
      os << endl;
    else
      os << " ]" << endl;
  }
}

void td_real::dump(const string &name, std::ostream &os) const {
  std::ios_base::fmtflags old_flags = os.flags();
  std::streamsize old_prec = os.precision(19);
  os << std::scientific;

  string::size_type len = name.length();
  if (len > 0) {
    os << name << " = ";
    len += 3;
  }
  os << "[ ";
  len += 2;
  os << setw(27) << x[0] << ", " << setw(26) << x[1] << "," << endl;
  for (string::size_type i = 0; i < len; i++) os << ' ';
  os << setw(27) << x[2] << "  ]" << endl;

  os.precision(old_prec);
  os.flags(old_flags);
}

/* Divisions */
/* triple-double / double-double */
td_real td_real::sloppy_div(const td_real &a, const dd_real &b) {
  double q0, q1, q2;
  td_real r;
  td_real td_b(b);

  q0 = a[0] / b._hi();
  r = a - q0 * td_b;

  q1 = r[0] / b._hi();
  r -= (q1 * td_b);

  q2 = r[0] / b._hi();

  renorm3(q0, q1, q2);
  return td_real(q0, q1, q2);
}

td_real td_real::accurate_div(const td_real &a, const dd_real &b) {
  double q0, q1, q2, q3;
  td_real r;
  td_real td_b(b);

  q0 = a[0] / b._hi();
  r = a - q0 * td_b;

  q1 = r[0] / b._hi();
  r -= (q1 * td_b);

  q2 = r[0] / b._hi();
  r -= (q2 * td_b);

  q3 = r[0] / b._hi();

  renorm3(q0, q1, q2, q3);
  return td_real(q0, q1, q2);
}

/* triple-double / triple-double */
td_real td_real::sloppy_div(const td_real &a, const td_real &b) {
  double q0, q1, q2;

  td_real r;

  q0 = a[0] / b[0];
  r = a - (b * q0);

  q1 = r[0] / b[0];
  r -= (b * q1);

  q2 = r[0] / b[0];

  renorm3(q0, q1, q2);

  return td_real(q0, q1, q2);
}

td_real td_real::accurate_div(const td_real &a, const td_real &b) {
  double q0, q1, q2;

  td_real r;

  q0 = a[0] / b[0];
  r = a - (b * q0);

  q1 = r[0] / b[0];
  r -= (b * q1);

  q2 = r[0] / b[0];

  r -= (b * q2);
  double q3 = r[0] / b[0];

  renorm3(q0, q1, q2, q3);

  return td_real(q0, q1, q2);
}

/* Long division driven by the branch-free triple-word FMA.  Same
   correction sequence as accurate_div, but each residual
   r <- r - q * b  is a single fused tw_fma instead of a multiply followed
   by a subtraction, which removes one renormalization per step. */
td_real td_real::fma_div(const td_real &a, const td_real &b) {
  double q0, q1, q2, q3;
  td_real r;

  q0 = a[0] / b[0];
  r = tw_fma_safe(b, -q0, a);       /* r = a - q0 * b */

  q1 = r[0] / b[0];
  r = tw_fma_safe(b, -q1, r);

  q2 = r[0] / b[0];
  r = tw_fma_safe(b, -q2, r);

  q3 = r[0] / b[0];

  renorm3(q0, q1, q2, q3);

  return td_real(q0, q1, q2);
}

QD_API td_real sqrt(const td_real &a) {
  /* Newton iteration:  r' = r + (0.5 - h * r^2) * r  with h = a/2,
     converging to 1/sqrt(a).  Starts with a double approximation.

     Two changes over the 0.0.2 version:
       * the first refinement runs in dd_real, which is already enough to
         carry ~32 digits, so only two triple-double steps are needed;
       * each step is a pair of branch-free tw_fma calls, so neither
         h * r^2 nor the following multiply-and-add is renormalized on
         its own.                                                        */

  if (a.is_zero())
    return 0.0;

  if (a.is_negative()) {
    td_real::error("(td_real::sqrt): Negative argument.");
    return td_real::_nan;
  }

  /* ~16 digits -> ~32 digits, in double-double. */
  dd_real ad(a[0], a[1]);
  dd_real hd = mul_pwr2(ad, 0.5);
  dd_real rd = (1.0 / std::sqrt(a[0]));
  rd = dw_fma(dw_fma(-hd, sqr(rd), dd_real(0.5)), rd, rd);

  td_real h = mul_pwr2(a, 0.5);
  td_real r(rd);

  /* ~32 digits -> full triple-double precision. */
  r = tw_fma(tw_fma(-h, sqr(r), td_real(0.5)), r, r);

  r *= a;
  return r;
}

/* Reference (pre-0.0.3) square root, kept so that the benchmark can
   measure what the branch-free FMA buys us. */
QD_API td_real sqrt_legacy(const td_real &a) {
  if (a.is_zero())
    return 0.0;

  if (a.is_negative()) {
    td_real::error("(td_real::sqrt_legacy): Negative argument.");
    return td_real::_nan;
  }

  td_real r = (1.0 / std::sqrt(a[0]));
  td_real h = mul_pwr2(a, 0.5);

  r += ((0.5 - h * sqr(r)) * r);
  r += ((0.5 - h * sqr(r)) * r);
  r += ((0.5 - h * sqr(r)) * r);

  r *= a;
  return r;
}

td_real nroot(const td_real &a, int n) {
  if (n <= 0) {
    td_real::error("(td_real::nroot): N must be positive.");
    return td_real::_nan;
  }

  if (n % 2 == 0 && a.is_negative()) {
    td_real::error("(td_real::nroot): Negative argument.");
    return td_real::_nan;
  }

  if (n == 1) {
    return a;
  }
  if (n == 2) {
    return sqrt(a);
  }
  if (a.is_zero()) {
    return td_real(0.0);
  }

  td_real r = abs(a);
  td_real x = std::exp(-std::log(r.x[0]) / n);

  /* Newton's iteration. */
  double dbl_n = static_cast<double>(n);
  x += x * (1.0 - r * npwr(x, n)) / dbl_n;
  x += x * (1.0 - r * npwr(x, n)) / dbl_n;
  x += x * (1.0 - r * npwr(x, n)) / dbl_n;
  if (a[0] < 0.0) {
    x = -x;
  }
  return 1.0 / x;
}

/* exp / expm1 / log / log10 / sin / cos / sincos / sinh / cosh have
   native triple-double implementations at the bottom of qd_real.cpp
   (they share that file's inv_fact and sin/cos tables by truncation).
   The remaining transcendentals below are still computed by promoting
   to qd_real, calling the qd_real implementation, and truncating back
   to td_real.                                                          */

td_real tan(const td_real &a) {
  return td_real(tan(to_qd_real(a)));
}

td_real asin(const td_real &a) {
  return td_real(asin(to_qd_real(a)));
}

td_real acos(const td_real &a) {
  return td_real(acos(to_qd_real(a)));
}

td_real atan(const td_real &a) {
  return td_real(atan(to_qd_real(a)));
}

td_real atan2(const td_real &y, const td_real &x) {
  return td_real(atan2(to_qd_real(y), to_qd_real(x)));
}

td_real tanh(const td_real &a) {
  return td_real(tanh(to_qd_real(a)));
}

void sincosh(const td_real &a, td_real &s, td_real &c) {
  qd_real qs, qc;
  sincosh(to_qd_real(a), qs, qc);
  s = td_real(qs);
  c = td_real(qc);
}

td_real asinh(const td_real &a) {
  return td_real(asinh(to_qd_real(a)));
}

td_real acosh(const td_real &a) {
  return td_real(acosh(to_qd_real(a)));
}

td_real atanh(const td_real &a) {
  return td_real(atanh(to_qd_real(a)));
}

td_real drem(const td_real &a, const td_real &b) {
  td_real n = nint(a/b);
  return (a - n * b);
}

td_real divrem(const td_real &a, const td_real &b, td_real &r) {
  td_real n = nint(a/b);
  r = a - n * b;
  return n;
}

QD_API td_real fmod(const td_real &a, const td_real &b) {
  td_real n = aint(a / b);
  return (a - b * n);
}

QD_API td_real tdrand() {
  static const double m_const = 4.6566128730773926e-10;  /* = 2^{-31} */
  double m = m_const;
  td_real r = 0.0;
  double d;

  /* Generate ~31 bits at a time, accumulate over 6 iterations
     for ~186 bits of randomness (more than td_real precision). */
  for (int i = 0; i < 6; i++, m *= m_const) {
    d = std::rand() * m;
    r += d;
  }

  return r;
}


/* polyeval(c, n, x) */
td_real polyeval(const td_real *c, int n, const td_real &x) {
  td_real r = c[n];

  for (int i = n-1; i >= 0; i--) {
    r *= x;
    r += c[i];
  }

  return r;
}

QD_API td_real polyroot(const td_real *c, int n,
    const td_real &x0, int max_iter, double thresh) {
  td_real x = x0;
  td_real f;
  td_real *d = new td_real[n];
  bool conv = false;
  int i;
  double max_c = std::abs(to_double(c[0]));
  double v;

  if (thresh == 0.0) thresh = td_real::_eps;

  for (i = 1; i <= n; i++) {
    v = std::abs(to_double(c[i]));
    if (v > max_c) max_c = v;
    d[i-1] = c[i] * static_cast<double>(i);
  }
  thresh *= max_c;

  for (i = 0; i < max_iter; i++) {
    f = polyeval(c, n, x);

    if (abs(f) < thresh) {
      conv = true;
      break;
    }
    x -= (f / polyeval(d, n-1, x));
  }
  delete [] d;

  if (!conv) {
    td_real::error("(td_real::polyroot): Failed to converge.");
    return td_real::_nan;
  }

  return x;
}

td_real td_real::debug_rand() {
  if (std::rand() % 2 == 0)
    return tdrand();

  int expn = 0;
  td_real a = 0.0;
  double d;
  for (int i = 0; i < 3; i++) {
    d = std::ldexp(std::rand() / static_cast<double>(RAND_MAX), -expn);
    a += d;
    expn = expn + 54 + std::rand() % 200;
  }
  return a;
}
