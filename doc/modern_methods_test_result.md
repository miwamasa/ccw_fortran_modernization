# 現代的手法のテスト結果 — VA10A vs Levenberg-Marquardt vs Broyden

## テスト環境

- **コンパイラ:** GNU Fortran (GCC) 13.3.0
- **ビルドシステム:** CMake 3.28.3
- **OS:** Linux 6.18.5

## テスト関数

| テスト | 関数 | 残差形式 | 解析解 |
|--------|------|---------|--------|
| Test 1 | Rosenbrock: f=(1-x)²+100(y-x²)² | r₁=1-x, r₂=10(y-x²) | (1, 1), f=0 |
| Test 2 | 二次関数: f=(x-3)²+(y+2)² | r₁=x-3, r₂=y+2 | (3, -2), f=0 |

## 結果比較

### Test 1: Rosenbrock関数（初期点: (-1, 1)）

| 手法 | x | y | f(x,y) | iexit | 結果 |
|------|---|---|--------|-------|------|
| VA10A (Original) | 1.000000 | 1.000000 | 2.98×10⁻¹⁷ | 1 (収束) | PASS |
| Levenberg-Marquardt | 1.000000 | 1.000000 | 4.55×10⁻²⁰ | 1 (収束) | PASS |
| Broyden | 1.000000 | 1.000000 | 4.93×10⁻³⁰ | 1 (収束) | PASS |

### Test 2: 二次関数（初期点: (0, 0)）

| 手法 | x | y | f(x,y) | iexit | 結果 |
|------|---|---|--------|-------|------|
| VA10A (Original) | 3.000000 | -2.000000 | 0.0 | 2 (勾配条件) | PASS |
| Levenberg-Marquardt | 3.000000 | -2.000000 | 2.31×10⁻²⁹ | 1 (収束) | PASS |
| Broyden | 3.000000 | -2.000000 | 0.0 | 1 (収束) | PASS |

### CTest 全体結果

```
Test project /home/user/ccw_fortran_modernization/build
    Start 1: minimize_test
1/3 Test #1: minimize_test ....................   Passed    0.00 sec
    Start 2: lm_test
2/3 Test #2: lm_test ..........................   Passed    0.00 sec
    Start 3: broyden_test
3/3 Test #3: broyden_test .....................   Passed    0.00 sec

100% tests passed, 0 tests failed out of 3
```

## 分析

### 精度の比較

Rosenbrock関数での最終関数値:
- VA10A: 2.98×10⁻¹⁷
- Levenberg-Marquardt: 4.55×10⁻²⁰ （VA10Aより約1000倍高精度）
- Broyden: 4.93×10⁻³⁰ （VA10Aより約10¹³倍高精度）

**Broyden法が最も高精度**。これは残差 r(x)=0 を直接解くため、二乗和を経由しないことが理由。

### 各手法の特性

| 特性 | VA10A | Levenberg-Marquardt | Broyden |
|------|-------|---------------------|---------|
| アルゴリズム種別 | 準Newton（一般最適化） | 非線形最小二乗特化 | 非線形方程式系ソルバー |
| 問題構造の活用 | なし（f(x)のみ使用） | 残差構造を活用（J^TJ ≈ H） | 残差を直接ゼロにする |
| Hessian近似 | Broyden族ランク2更新 | J^TJ + λI（Gauss-Newton + damping） | Jacobian逆行列のランク1更新 |
| 勾配計算 | 数値微分（前進/中心差分） | 数値Jacobian（前進差分） | 初期Jacobianのみ数値計算 |
| GOTO文 | 約30箇所 | なし（構造化制御フロー） | なし（構造化制御フロー） |
| コードの可読性 | 非常に低い | 高い | 高い |
| Fortran規格 | FORTRAN 77 | Fortran 2008+ | Fortran 2008+ |

### DMFTバスフィッティングへの適用

1. **Levenberg-Marquardt法**: DMFTのバスフィッティング問題は本質的に非線形最小二乗問題であるため、最も自然な選択。残差 |G₀(iωₙ) - G₀ᴬⁿᵈ(iωₙ)| を個別に扱えるため、Jacobian J^T J による Hessian 近似が問題構造に適合する。

2. **Broyden法**: DMFTの自己無撞着条件 G₀ᴬⁿᵈ = G₀ を不動点方程式として直接解く。特に松原周波数の数とパラメータの数が等しい場合（正方系）に最適。Johnson (1988) の混合法として DMFT コミュニティで広く使用されている。

## 結論

3つの手法はすべて同一の最適化問題を正しく解くことが確認された。現代的手法（LM, Broyden）はVA10Aと同等以上の精度を持ち、かつ構造化されたコード（GOTO文なし、モジュール化、明示的型宣言）で実装されている。DMFTバスフィッティングのモダナイゼーションにおいて、いずれも有力な代替手法である。
