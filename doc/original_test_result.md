# 初期テスト結果 — minimize サブルーチン単体テスト

## テスト環境

- **コンパイラ:** GNU Fortran (GCC) 13.3.0
- **ビルドシステム:** CMake 3.28.3
- **OS:** Linux 6.18.5

## ビルド手順

```bash
mkdir build && cd build
cmake ..
make
```

ビルドは成功。Fortran 2018で削除された機能（算術IF文、共有DOラベル等）に関する警告が出るが、コンパイル・リンクは正常完了。

## テスト実行結果

```
================================================
 minimize subroutine test suite
 ================================================

 Test 1: Rosenbrock function
   f(x,y) = (1-x)^2 + 100*(y-x^2)^2
   Expected minimum: f(1, 1) = 0

  Result: x =     1.000000, y =     1.000000
  f(x,y) =   0.2980800E-16
  iexit  =   1

   >>> PASS: Converged to correct minimum

 Test 2: Quadratic function
   f(x,y) = (x-3)^2 + (y+2)^2
   Expected minimum: f(3, -2) = 0

  Result: x =     3.000000, y =    -2.000000
  f(x,y) =   0.0000000E+00
  iexit  =   2

   >>> PASS: Converged to correct minimum

 ================================================
 Results:  2 PASS,  0 FAIL
 ================================================
```

## CTest結果

```
Test project /home/user/ccw_fortran_modernization/build
    Start 1: minimize_test
1/1 Test #1: minimize_test ....................   Passed    0.00 sec

100% tests passed, 0 tests failed out of 1

Total Test time (real) =   0.00 sec
```

## 結果分析

### テスト1: Rosenbrock関数

| 項目 | 値 |
|------|-----|
| 初期点 | (-1.0, 1.0) |
| 収束先 | (1.000000, 1.000000) |
| 解析解 | (1.0, 1.0) |
| 最終関数値 | 2.98 × 10⁻¹⁷ (≈ 0) |
| iexit | 1 (収束) |

Rosenbrock関数は急な谷構造を持つ困難な最適化問題だが、minimizeは正確に最小点に収束した。iexit=1は正常収束を示す。

### テスト2: 二次関数

| 項目 | 値 |
|------|-----|
| 初期点 | (0.0, 0.0) |
| 収束先 | (3.000000, -2.000000) |
| 解析解 | (3.0, -2.0) |
| 最終関数値 | 0.0 (完全一致) |
| iexit | 2 (勾配条件で停止) |

単純な二次関数では厳密解に到達。iexit=2は勾配の内積が非負になった（完全な最小点に到達した）ことを示す。

## 結論

オリジナルのFortranコードのminimizeサブルーチンは、両方のテスト関数に対して正しく動作することが確認された。この結果は、モダナイゼーション後の等価性検証のベースラインとして使用する。
