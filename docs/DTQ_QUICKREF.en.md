# dtq Quick Reference (C++ / C Usage Guide)

A short reference for using `dtq` (a CPU-side multi-precision floating-point
library derived from QD 2.13) from your own project.

Target version: **dtq-0.0.2**

For the GPU (CUDA) counterpart, see the separate `gdtq` package.

---

## 1. Types

`dtq` represents multi-precision values as **unevaluated sums** of IEEE-754
values. Two families are provided: a double-based one (the usual choice) and
a float-based one (added in 0.0.2, useful as a SIMD- or GPU-style lightweight
container).

### 1.1 Double-based (recommended)

| Type | Layout | Precision | Notes |
|---|---|---|---|
| `dd_real` | `double[2]` | ~32 digits (2×53 bit, ~106 bit) | DD; from QD |
| `td_real` | `double[3]` | ~47 digits (3×53 bit, ~156 bit) | TD; new in this package |
| `qd_real` | `double[4]` | ~63 digits (4×53 bit, ~212 bit) | QD |

Headers: `<qd/dd_real.h>`, `<qd/td_real.h>`, `<qd/qd_real.h>`

### 1.2 Float-based (added in dtq-0.0.2)

| Type | Layout | Precision | Notes |
|---|---|---|---|
| `ds_real` | `float[2]` | ~14 digits | DS; ~double precision in two floats |
| `ts_real` | `float[3]` | ~21 digits | TS |
| `qs_real` | `float[4]` | ~28 digits | QS |

Headers: `<qd/ds_real.h>`, `<qd/ts_real.h>`, `<qd/qs_real.h>`

---

## 2. Build and install

```sh
./configure
make
make check         # runs qd_test, c_test, pslq_test
sudo make install
```

Useful `configure` options:

| Option | Meaning |
|---|---|
| `--prefix=DIR` | install prefix (default `/usr/local`) |
| `--enable-fma=auto` | use FMA when available (default) |
| `--enable-fma=no` | disable FMA |
| `--enable-debug` | `-g` and extra debug checks |
| `--disable-shared` | build only the static library |

### Important: floating-point flags

For GCC/Clang, **`-ffp-contract=off` is added automatically**. Do not strip
it: the error-free transformations used by dd/td/qd (Two-Sum / Two-Prod / ...)
break if the compiler is allowed to silently fuse `a*b + c` into an FMA.
With FMA fusion in place, exp/log lose ~12 digits.

On x86_64, `-mfma` is also added automatically so that the explicit
`__builtin_fma` calls inside the library compile to a hardware FMA3
instruction. AArch64 / IA-64 / PowerPC have FMA in the base ISA and need no
extra flag.

### Regenerating `configure`

```sh
autoreconf -fi
```

---

## 3. Using the library from C++

### 3.1 Minimal example

```cpp
#include <iostream>
#include <iomanip>
#include <qd/dd_real.h>
#include <qd/td_real.h>
#include <qd/qd_real.h>
#include <qd/fpu.h>

int main() {
    unsigned int cw;
    fpu_fix_start(&cw);                          // pin x86 FPU to 64-bit

    dd_real a = dd_real::_pi;                    // a constant
    td_real b = sqrt(td_real(2.0));              // a math function
    qd_real c = exp(qd_real(1.0));
    qd_real d = qd_real("1.234567890123456789012345678901234567890");

    std::cout << std::setprecision(dd_real::_ndigits) << "pi = " << a << "\n";
    std::cout << std::setprecision(td_real::_ndigits) << "sqrt(2) = " << b << "\n";
    std::cout << std::setprecision(qd_real::_ndigits) << "e = " << c << "\n";

    fpu_fix_end(&cw);
    return 0;
}
```

### 3.2 Compiling

The recommended way is to use `qd-config`, which supplies the right flags
(including `-ffp-contract=off`):

```sh
g++ `qd-config --cxxflags` example.cpp `qd-config --libs` -o example
```

By hand:

```sh
g++ -O2 -ffp-contract=off -I/usr/local/include \
    example.cpp -L/usr/local/lib -lqd -lm -o example
```

### 3.3 Available operations

| Category | Examples |
|---|---|
| Arithmetic | `+ - * /`, `+= -= *= /=` (`X_real op X_real`, `X_real op double`, `double op X_real`) |
| Comparison | `== != < <= > >=` |
| Unary | `-a`, `abs(a)`, `sqr(a)`, `sqrt(a)`, `nroot(a, n)`, `npwr(a, n)` |
| Rounding | `nint(a)`, `floor(a)`, `ceil(a)`, `aint(a)` |
| Exp / log | `exp(a)`, `log(a)`, `log10(a)` |
| Trig | `sin`, `cos`, `tan`, `asin`, `acos`, `atan`, `atan2(y, x)`, `sincos(a, s, c)` |
| Hyperbolic | `sinh`, `cosh`, `tanh`, `asinh`, `acosh`, `atanh` |
| I/O | `std::cout << a`, `std::cin >> a`, `a.to_string(precision)`, `X_real("123.456")` |
| Predicates | `a.is_zero()`, `a.is_one()`, `a.is_positive()`, `a.is_negative()`, `a.isnan()`, `a.isinf()`, `a.isfinite()` |

### 3.4 Constants (same names across dd/td/qd)

```
X_real::_pi   X_real::_2pi   X_real::_pi2   X_real::_pi4   X_real::_3pi4
X_real::_e    X_real::_log2  X_real::_log10
X_real::_eps  X_real::_max   X_real::_safe_max  X_real::_min_normalized
X_real::_nan  X_real::_inf
X_real::_ndigits         // number of decimal digits to use with setprecision
```

### 3.5 The float-based types follow the same shape

```cpp
#include <qd/ds_real.h>
ds_real x = ds_real("3.14159265358979");
ds_real y = sqrt(x) * x;
std::cout << std::setprecision(ds_real::_ndigits) << y << "\n";
```

For ds/ts/qs, `_eps` is a `float`, and constants such as `_pi` are of the
respective `ds_real` / `ts_real` / `qs_real` type.

---

## 4. Using the library from C

C wrappers are provided for the **double-based** types (`dd`, `td`, `qd`).
The float-based types (`ds`, `ts`, `qs`) **do not** have a C wrapper.

### 4.1 Headers and representation

| Header | Storage |
|---|---|
| `<qd/c_dd.h>` | `double[2]` |
| `<qd/c_td.h>` | `double[3]` |
| `<qd/c_qd.h>` | `double[4]` |

### 4.2 Minimal example

```c
#include <stdio.h>
#include <qd/c_qd.h>
#include <qd/fpu.h>

int main(void) {
    fpu_fix_start(NULL);

    double a[4], b[4], c[4];

    c_qd_copy_d(2.0, a);          /* a = 2.0           */
    c_qd_sqrt(a, b);              /* b = sqrt(2)       */
    c_qd_pi(c);                   /* c = pi            */
    c_qd_mul(b, c, a);            /* a = sqrt(2) * pi  */
    c_qd_write(a);                /* print to stdout   */

    return 0;
}
```

The most thorough example is `tests/c_test.c`, which uses the Salamin–Brent
formula to compute π in dd / td / qd.

### 4.3 Naming convention

```
c_<prec>_<op>[ _<argtypes>]
```

`<prec>` is `dd` / `td` / `qd`. `<op>` is one of `add`, `sub`, `mul`, `div`,
`sqr`, `sqrt`, `exp`, `log`, `sin`, `cos`, `tan`, `pi`, `copy`, `comp` (compare),
…. Compound assignments take the form `c_<prec>_self<op>` (for example,
`c_qd_selfmul_d(2.0, p)` is `p *= 2.0`).

Mixed-precision helpers also exist (e.g. `c_qd_add_dd_qd` adds a DD and a QD
into a QD). See the headers for the full list.

### 4.4 Compiling and linking

```sh
gcc -O2 -ffp-contract=off -I/usr/local/include \
    mycode.c -L/usr/local/lib -lqd -lstdc++ -lm -o mycode
```

The C wrappers are implemented in C++, so linking needs `libstdc++`
(or `libc++`). Either pass `-lstdc++` to the C linker, or simply link with
`g++`.

---

## 5. Reference programs

| Program | Purpose |
|---|---|
| `tests/qd_test` | sanity / accuracy checks for dd/td/qd math functions |
| `tests/c_test` | exercises the C wrappers (Salamin–Brent for π) |
| `tests/pslq_test` | PSLQ integer-relation search at all precisions |
| `tests/qd_timer` | simple benchmark (built via `make demo`) |
| `tests/quadt_test` | tanh-sinh quadrature (built via `make demo`) |
| `tests/huge` | huge-number formatting demo (`make demo`) |

Examples:

```sh
cd tests
./qd_test -v        # verbose output for every precision
./qd_test -td       # only triple-double
./pslq_test -all -v
```

---

## 6. Common pitfalls

| Symptom | Cause / fix |
|---|---|
| exp/log lose ~10¹² eps on ARM64 / modern x86_64 | Compiler fused mul+add into FMA via `-ffp-contract=fast`. Always pass `-ffp-contract=off` (auto-added under GCC/Clang in this package, but easy to lose if you set `CXXFLAGS=...` yourself) |
| Results jitter on x86 | x87 80-bit extended precision is on. Call `fpu_fix_start(&cw)` at the top of `main`. |
| Undefined symbols when linking C code | Link with `g++`, or pass `-lstdc++ -lm` after `-lqd` |
| `qd_real("0.1")` is slow in a hot loop | Decimal→binary conversion runs every time. Hoist the constant out of the loop. |
| Output looks shorter than `_ndigits` would suggest | `std::setprecision(_ndigits)` must be in effect on every output, or set persistently via `std::cout.precision(...)` |

---

## 7. Files worth reading

- `tests/qd_test.cpp` — assorted C++ usage examples
- `tests/c_test.c` — C-wrapper usage (Salamin–Brent)
- `tests/pslq_test.cpp` — a non-trivial algorithm implementation
- `include/qd/dd_real.h` / `td_real.h` / `qd_real.h` — full C++ API
- `include/qd/c_dd.h` / `c_td.h` / `c_qd.h` — full C wrapper API
- `docs/qd.pdf` — original QD library theory paper (LaTeX build)
