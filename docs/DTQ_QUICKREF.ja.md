# dtq クイックリファレンス（C++ / C 利用ガイド）

CPU 上で多倍長浮動小数点演算を行うライブラリ `dtq`（QD 2.13 派生）を、
自分のプロジェクトから利用するときの簡易リファレンスです。

対象バージョン: **dtq-0.0.3**

GPU（CUDA）版は別パッケージ `gdtq` を参照してください。

---

## 1. 提供される型

`dtq` は IEEE-754 値の **unevaluated sum** で多倍長を表現します。
double を使う系統と float を使う系統の 2 系統があります。

### 1.1 double 系（推奨：通常はこちら）

| 型 | 構造 | 概略精度 | 主な用途 |
|---|---|---|---|
| `dd_real` | `double[2]` | ~32 桁 (2×53bit, ~106bit) | DD：QD 由来、最も使用例が多い |
| `td_real` | `double[3]` | ~47 桁 (3×53bit, ~156bit) | TD：本パッケージで追加 |
| `qd_real` | `double[4]` | ~63 桁 (4×53bit, ~212bit) | QD |

ヘッダ: `<qd/dd_real.h>` `<qd/td_real.h>` `<qd/qd_real.h>`

### 1.2 float 系（dtq-0.0.2 で追加。SIMD/GPU 寄りの軽量型）

| 型 | 構造 | 概略精度 | 用途 |
|---|---|---|---|
| `ds_real` | `float[2]` | ~14 桁 | DS：double 相当を float 2本で |
| `ts_real` | `float[3]` | ~21 桁 | TS |
| `qs_real` | `float[4]` | ~28 桁 | QS |

ヘッダ: `<qd/ds_real.h>` `<qd/ts_real.h>` `<qd/qs_real.h>`

**有効値域に注意**: float 系の値 x はリム（構成要素）が概ね
|x|, |x|·2⁻²⁴, |x|·2⁻⁴⁸, |x|·2⁻⁷² の位置に並ぶため、最下位リムが
float の正規化数下限 2⁻¹²⁶ を割ると精度が落ち始めます。目安として

| 型 | フル精度が保てる範囲 |
|---|---|
| `ds_real` | \|x\| ≥ 2⁻¹⁰² (約 2.0e-31) |
| `ts_real` | \|x\| ≥ 2⁻⁷⁸ (約 3.3e-24) |
| `qs_real` | \|x\| ≥ 2⁻⁵⁴ (約 5.5e-17) |

これはアルゴリズムではなく形式そのものの表現限界です（例：商がしきい値を
下回る除算では誤差が数十 eps 以上になり得ます）。double 系では対応する
しきい値が 2⁻⁸⁶⁶ 以下なので実用上問題になりません。

---

## 2. ビルドとインストール

```sh
./configure
make
make check         # qd_test, c_test, pslq_test, fma_test 実行
sudo make install
```

主な `configure` オプション:

| オプション | 意味 |
|---|---|
| `--prefix=DIR` | インストール先（既定 `/usr/local`） |
| `--enable-fma=auto` | FMA を可能なら有効化（既定） |
| `--enable-fma=no` | FMA を無効化 |
| `--enable-debug` | `-g` とデバッグチェックを有効化 |
| `--disable-shared` | 静的ライブラリのみ |

### 重要：浮動小数点フラグ

GCC/Clang では **自動で `-ffp-contract=off` が付与** されます。これは絶対に
外さないでください。Two-Sum / Two-Prod などの誤差なし変換は、コンパイラが
`a*b + c` を FMA に融合すると壊れます（exp/log で約 12 桁失います）。

x86_64 では `-mfma` も自動付与されます（ライブラリ内で明示的に呼んでいる
`__builtin_fma` をハードウェア FMA3 命令に落とすため）。AArch64 / IA-64 /
PowerPC は基本 ISA に FMA があるので追加フラグ不要です。

### `configure` を再生成するとき

```sh
autoreconf -fi
```

---

## 3. C++ から使う

### 3.1 最小例

```cpp
#include <iostream>
#include <iomanip>
#include <qd/dd_real.h>
#include <qd/td_real.h>
#include <qd/qd_real.h>
#include <qd/fpu.h>

int main() {
    unsigned int cw;
    fpu_fix_start(&cw);                          // x86 の FPU を 64bit に固定

    dd_real a = dd_real::_pi;                    // 定数
    td_real b = sqrt(td_real(2.0));              // 数学関数
    qd_real c = exp(qd_real(1.0));
    qd_real d = qd_real("1.234567890123456789012345678901234567890"); // 文字列から

    std::cout << std::setprecision(dd_real::_ndigits) << "pi = " << a << "\n";
    std::cout << std::setprecision(td_real::_ndigits) << "sqrt(2) = " << b << "\n";
    std::cout << std::setprecision(qd_real::_ndigits) << "e = " << c << "\n";

    fpu_fix_end(&cw);
    return 0;
}
```

### 3.2 ビルド

`qd-config` を使うのが推奨（`-ffp-contract=off` などのフラグを正しく入れて
くれる）:

```sh
g++ `qd-config --cxxflags` example.cpp `qd-config --libs` -o example
```

`pkg-config` 風に手書きする場合:

```sh
g++ -O2 -ffp-contract=off -I/usr/local/include \
    example.cpp -L/usr/local/lib -lqd -lm -o example
```

### 3.3 主な操作

| カテゴリ | 例 |
|---|---|
| 算術 | `+ - * /`、`+= -= *= /=`（`X_real op X_real`、`X_real op double`、`double op X_real`） |
| 比較 | `== != < <= > >=` |
| 単項 | `-a`、`abs(a)`、`sqr(a)`、`sqrt(a)`、`nroot(a, n)`、`npwr(a, n)` |
| 積和（FMA） | `fma(a, b, c)`（= `a*b + c`）、型別の名前は `dw_fma`（dd/ds）、`tw_fma`（td/ts）、`qw_fma`（qd/qs） |
| 丸め | `nint(a)`、`floor(a)`、`ceil(a)`、`aint(a)` |
| 指数・対数 | `exp(a)`、`log(a)`、`log10(a)` |
| 三角関数 | `sin(a)`、`cos(a)`、`tan(a)`、`asin(a)`、`acos(a)`、`atan(a)`、`atan2(y, x)`、`sincos(a, s, c)` |
| 双曲線 | `sinh(a)`、`cosh(a)`、`tanh(a)`、`asinh(a)`、`acosh(a)`、`atanh(a)` |
| 入出力 | `std::cout << a`、`std::cin >> a`、`a.to_string(precision)`、`X_real("123.456")` |
| 述語 | `a.is_zero()`、`a.is_one()`、`a.is_positive()`、`a.is_negative()`、`a.isnan()`、`a.isinf()`、`a.isfinite()` |

### 3.3b 分岐なし積和演算（fused multiply-add, 0.0.3 で追加）

`a * b + c` を 1 個の融合演算として計算するルーチンを全精度クラスに
用意しています。

```cpp
dd_real z = dw_fma(a, b, c);   // double-word （dd_real, ds_real）
td_real z = tw_fma(a, b, c);   // triple-word （td_real, ts_real）
qd_real z = qw_fma(a, b, c);   // quad-word   （qd_real, qs_real）

qd_real z = fma(a, b, c);      // 総称名。6 型すべてで使えます
qd_real z = qw_fma(a, 2.5, c); // 乗数は素の double / float でも可
```

積 `a*b` の各項と加数 `c` の各語を 1 本の直線的（分岐なし）な加算
ネットワークにまとめて累算するため、`a*b` を単独で正規化しません。
`a * b + c` と書くより演算数が少なく、精度が落ちることもありません。

除算と平方根はこのルーチンの上に構築されています。`operator/` は
`T::fma_div` を呼びます（`-DQD_NO_FMA_DIV` を付けると 0.0.2 の
`sloppy_div` / `accurate_div` に戻ります）。`sqrt` は融合積和による
Newton 反復を使います。0.0.2 の実装は `sqrt_legacy(a)` として残して
あります。

手元の環境での新旧比較は `tests/` で `make bench` を実行してください。

### 3.4 主な定数（dd/td/qd 共通の命名）

```
X_real::_pi   X_real::_2pi   X_real::_pi2   X_real::_pi4   X_real::_3pi4
X_real::_e    X_real::_log2  X_real::_log10
X_real::_eps  X_real::_max   X_real::_safe_max  X_real::_min_normalized
X_real::_nan  X_real::_inf
X_real::_ndigits         // 小数点以下桁数（ストリーム出力時の precision に）
```

### 3.5 float 系（ds/ts/qs）の使い方も同形

```cpp
#include <qd/ds_real.h>
ds_real x = ds_real("3.14159265358979");
ds_real y = sqrt(x) * x;
std::cout << std::setprecision(ds_real::_ndigits) << y << "\n";
```

ds/ts/qs は `_eps` が `float` で、`_pi` などが各々 `ds_real` / `ts_real` /
`qs_real` 型で提供されます。

---

## 4. C から使う

dd / td / qd には C ラッパが用意されています（float 系の ds/ts/qs には
C ラッパは **ありません**）。

### 4.1 ヘッダと表現

| ヘッダ | 型の表現 |
|---|---|
| `<qd/c_dd.h>` | `double[2]` |
| `<qd/c_td.h>` | `double[3]` |
| `<qd/c_qd.h>` | `double[4]` |

### 4.2 最小例

```c
#include <stdio.h>
#include <qd/c_qd.h>
#include <qd/fpu.h>

int main(void) {
    fpu_fix_start(NULL);

    double a[4], b[4], c[4];

    c_qd_copy_d(2.0, a);          /* a = 2.0 */
    c_qd_sqrt(a, b);              /* b = sqrt(2)        */
    c_qd_pi(c);                   /* c = pi             */
    c_qd_mul(b, c, a);            /* a = sqrt(2) * pi   */
    c_qd_write(a);                /* 標準出力に表示     */

    return 0;
}
```

`tests/c_test.c` にある Salamin–Brent 公式の例が一番充実した参考例です。

### 4.3 関数命名規則

```
c_<prec>_<op>[ _<argtypes>]
```

`<prec>` は `dd` / `td` / `qd`、`<op>` は `add`、`sub`、`mul`、`div`、`sqr`、`sqrt`、`fma`、
`exp`、`log`、`sin`、`cos`、`tan`、`pi`、`copy`、`comp`（比較）など。
代入演算は `c_<prec>_self<op>` 形式（例: `c_qd_selfmul_d(2.0, p)` で `p *= 2.0`）。

混合演算用に `c_qd_add_dd_qd`（DD と QD を加算して QD を返す）等の補助も
あります。詳細は `<qd/c_qd.h>` 等のヘッダを参照してください。

### 4.4 ビルドとリンク

```sh
gcc -O2 -ffp-contract=off -I/usr/local/include \
    mycode.c -L/usr/local/lib -lqd -lstdc++ -lm -o mycode
```

C ラッパの実装は **C++ で書かれている** ため、リンクには `libstdc++`
（または `libc++`）が必要です。`gcc` リンカに `-lstdc++` を渡すか、`g++`
でリンクしてください。

---

## 5. テストプログラム（参考実装）

| プログラム | 内容 |
|---|---|
| `tests/qd_test` | dd/td/qd の数学関数の精度確認 |
| `tests/c_test` | C ラッパの動作確認（Salamin–Brent で π を計算） |
| `tests/pslq_test` | PSLQ 整数関係探索アルゴリズム |
| `tests/fma_test` | `dw_fma`/`tw_fma`/`qw_fma` と、それを使う除算・平方根の精度確認（`make check` で実行） |
| `tests/fma_bench` | 上記の新旧速度比較（`make bench`） |
| `tests/qd_timer` | 簡易ベンチマーク（`make demo` で構築） |
| `tests/quadt_test` | tanh-sinh 数値積分（`make demo`） |
| `tests/huge` | 大きな値の出力デモ（`make demo`） |

実行例:

```sh
cd tests
./qd_test -v        # 全精度の verbose 出力
./qd_test -td       # td_real のみ
./pslq_test -all -v
```

---

## 6. よくある落とし穴

| 症状 | 原因と対処 |
|---|---|
| ARM64 や AVX2 マシンで exp/log の精度が ~10¹² eps 落ちる | コンパイラが `-ffp-contract=fast` で FMA 融合した。`-ffp-contract=off` を必ず付ける（GCC/Clang なら自動だが、ユーザコード側でも忘れずに） |
| x86 で 80bit 拡張精度のせいで結果がぶれる | `fpu_fix_start(&cw)` を main 先頭で呼ぶ |
| C コードのリンクで未定義シンボル | `g++` でリンク、または `gcc -lqd -lstdc++ -lm` |
| 文字列構築が遅い | `qd_real("0.1")` は内部で十進→二進変換が走る。ループの外に出す |
| `_ndigits` で精度を指定しているのに出力が短い | `std::setprecision()` を毎回呼ぶ。`std::cout.precision()` の永続設定でも可 |

---

## 7. 参考ファイル

- `tests/qd_test.cpp` … C++ 側の各種使用例
- `tests/c_test.c` … C ラッパの使用例（Salamin–Brent）
- `tests/pslq_test.cpp` … 実用的なアルゴリズム実装例
- `include/qd/dd_real.h` `td_real.h` `qd_real.h` … C++ クラスの完全な API
- `include/qd/c_dd.h` `c_td.h` `c_qd.h` … C ラッパの完全な API
- `docs/qd.pdf` … 元の QD ライブラリの理論ペーパー（要 LaTeX ビルド）
