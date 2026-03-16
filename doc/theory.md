# 理論構築 — minimizeサブルーチンのモダナイゼーション基盤

## 1. 序論：なぜ単純なリファクタリングでは不十分か

lisalanc.fにおける`minimize`サブルーチン（VA10A）のモダナイゼーションを考えるとき、GOTO文の排除やFortran90構文への書き換えといった表層的なリファクタリングだけでは、このコードが持つ本質的な問題に対処できない。

VA10Aは1970年代のHarwellライブラリに由来する「ブラックボックス」であり、原作者も不明とされる。このルーチンが「なぜ動くのか」を理論的に再構築し、その数学的基盤の上に現代的な等価プログラムを設計することが、真のモダナイゼーションである。

## 2. VA10Aの数学的再構築

### 2.1 準Newton法の一般理論

無制約最適化問題 $\min_{x \in \mathbb{R}^n} f(x)$ に対して、準Newton法は以下の反復を行う：

$$x_{k+1} = x_k + \alpha_k d_k$$

ここで探索方向 $d_k = -B_k^{-1} \nabla f(x_k)$ であり、$B_k$ はHessian行列 $\nabla^2 f(x_k)$ の近似である。

VA10Aは $B_k$ を直接保持するのではなく、そのCholesky分解 $B_k = L_k D_k L_k^T$ を保持し、パック対称行列形式で格納する。これにより：
- 探索方向の計算が前進・後退代入で O(n²) で実行可能
- Hessian更新がCholesky因子のランク1更新として効率的に実行可能

### 2.2 Broyden族の更新公式

VA10Aのコア（行1099-1154）は、Broyden族のランク2更新を実装している。

一般のBroyden族更新は：

$$B_{k+1} = B_k + \frac{y_k y_k^T}{y_k^T s_k} - \frac{B_k s_k s_k^T B_k}{s_k^T B_k s_k} + \phi_k (s_k^T B_k s_k) v_k v_k^T$$

ここで：
- $s_k = x_{k+1} - x_k = \alpha_k d_k$（ステップ）
- $y_k = \nabla f(x_{k+1}) - \nabla f(x_k)$（勾配の変化）
- $\phi_k \in [0, 1]$: Broydenパラメータ

VA10Aのコードを分析すると：
- `link=1`（行1099-1113）: 第1ランク更新。$\phi_k$ の値に応じて BFGS型（$\phi=0$）または DFP型（$\phi=1$）の中間の更新
- `link=2`（行1114-1123）: 第2ランク更新。Broyden族の完全なランク2更新を構成
- 行1124-1154: Cholesky因子の直接更新（行列の再分解を回避）

### 2.3 線探索の理論

VA10Aの線探索（行1026-1082）は以下の戦略を組み合わせている：

1. **ステップ倍増（Doubling）:** 初期ステップで改善が見られる場合、ステップ幅を2倍にして探索（行1042-1054）
2. **安全停止条件:** `7f₁ + 5f₂ > 12f` による補間品質チェック（行1050-1051）。これは放物線補間の信頼性を判定する条件
3. **二次補間（Quadratic interpolation）:** 3点 $(0, f)$, $(\alpha, f_2)$, $(2\alpha, f_1)$ から放物線の最小点を推定（行1072-1074）

この線探索は厳密なWolfe条件を満たすものではないが、実用上十分なステップ幅削減を保証する。

## 3. 物理学的文脈：DMFT自己無撞着条件の変分的解釈

### 3.1 問題の定式化

lisalanc.fにおいて`minimize`が解く問題は、Anderson不純物模型のパラメータ $\{\epsilon_k, V_k\}$ を最適化して、Hubbardモデルから計算した裸のGreen関数 $\mathcal{G}_0(i\omega_n)$ とAnderson模型の $\mathcal{G}_0^{\text{And}}(i\omega_n)$ の不一致を最小化することである：

$$\min_{\{\epsilon_k, V_k\}} \mathcal{F}[\{\epsilon_k, V_k\}] = \frac{1}{N_\omega} \sum_{n=0}^{N_\omega} \left| \mathcal{G}_0(i\omega_n) - \mathcal{G}_0^{\text{And}}(i\omega_n; \{\epsilon_k, V_k\}) \right|$$

ここで Anderson模型のGreen関数は：

$$\mathcal{G}_0^{\text{And}}(i\omega_n) = \left[ i\omega_n - \epsilon_1 - \frac{U}{2} - \sum_{k=1}^{N_s-1} \frac{V_k^2}{i\omega_n - \epsilon_{k+1}} \right]^{-1}$$

### 3.2 関数空間上の変分問題としての再解釈

この最適化問題は、より高い視座から**関数空間上の変分問題**として再定式化できる。

**Hilbert空間 $\mathcal{H}$** を松原周波数上のGreen関数の空間とし、内積を：

$$\langle G_1, G_2 \rangle = \frac{1}{N_\omega} \sum_{n=0}^{N_\omega} G_1^*(i\omega_n) G_2(i\omega_n)$$

と定義する。このとき、目的汎関数は $\mathcal{H}$ 上のノルムで表現される：

$$\mathcal{F} = \| \mathcal{G}_0 - \mathcal{G}_0^{\text{And}} \|_{\mathcal{H}}$$

### 3.3 Galerkin法との対応

Anderson模型のGreen関数 $\mathcal{G}_0^{\text{And}}$ は、有限次元パラメータ空間 $\mathcal{P} = \{(\epsilon_k, V_k)\}_{k=1}^{N_s}$ によって張られるGreen関数の部分多様体 $\mathcal{M} \subset \mathcal{H}$ 上の元である。

自己無撞着条件は：

$$\text{Find } g^* \in \mathcal{M} \text{ such that } \| \mathcal{G}_0 - g^* \| = \min_{g \in \mathcal{M}} \| \mathcal{G}_0 - g \|$$

これは**Galerkin法**の変分定式化と正確に対応する：

- **試行関数空間:** Anderson模型パラメータで張られるGreen関数の集合 $\mathcal{M}$
- **残差直交条件:** 最適点 $g^*$ において、残差 $r = \mathcal{G}_0 - g^*$ は $\mathcal{M}$ の接空間に直交する

$$\langle r, \frac{\partial g^*}{\partial \epsilon_k} \rangle = 0, \quad \langle r, \frac{\partial g^*}{\partial V_k} \rangle = 0$$

この条件は、まさに目的関数の勾配をゼロにする条件に一致する。すなわち、VA10Aが数値的に解いている問題は、**Green関数空間上のGalerkin射影問題**に他ならない。

### 3.4 Ritz変分原理との関係

さらに深い視座として、DMFTの自己無撞着条件自体がLuttinger-Ward汎関数：

$$\Omega[\mathcal{G}] = -\text{Tr}\ln(-\mathcal{G}^{-1}) - \text{Tr}(\Sigma \mathcal{G}) + \Phi[\mathcal{G}]$$

の停留点条件から導かれる。Anderson模型パラメータ空間への制限は、この汎関数のRitz変分近似に対応する。

## 4. モダナイゼーションの理論的基盤

### 4.1 現代的最適化フレームワーク

VA10Aの理論的基盤を理解した上で、以下の現代的アルゴリズムが等価な代替として考えられる：

#### L-BFGS法（Limited-memory BFGS）

VA10Aが $O(n^2)$ のHessian記憶を必要とするのに対し、L-BFGSは最近 $m$ ステップの $\{s_k, y_k\}$ ペアのみを保持し、$O(mn)$ の記憶量で準Newton方向を近似する。

**利点:**
- 大規模問題（$n > 1000$）に対してメモリ効率的
- Nocedal (1980) の二重ループ再帰で効率的に実装可能
- VA10Aと同じBFGS族の更新に基づく

**等価性:** $m = n$ の場合、L-BFGSは完全BFGSと等価であり、VA10Aの挙動を再現する。

#### 信頼領域法（Trust Region Method）

$$\min_{d} m_k(d) = f_k + \nabla f_k^T d + \frac{1}{2} d^T B_k d \quad \text{s.t.} \quad \|d\| \leq \Delta_k$$

**利点:**
- 線探索が不要（信頼領域半径の調整で代替）
- 大域的収束が理論的に保証される
- Rosenbrock関数のような急な谷でもロバスト

#### 自動微分（Automatic Differentiation）

VA10Aの最大の弱点の一つは数値微分への依存である。現代のFortranコンパイラとライブラリは自動微分をサポートしており：

- **前進モードAD:** 各変数に対するタンジェントベクトルを伝播（$O(n)$ 回の関数評価相当）
- **逆向モードAD:** 1回の逆方向パスで全勾配を計算（$O(1)$ 回の関数評価相当）

数値微分の精度問題（前進差分の $O(h)$ 誤差、中心差分の $O(h^2)$ 誤差）を完全に排除できる。

### 4.2 Hilbert空間での最急降下法としての再定式化

Green関数空間 $\mathcal{H}$ 上での問題を直接扱うフレームワーク：

1. **Fréchet微分:** 目的汎関数のGreen関数に関する関数微分を定義
2. **Riemannian勾配降下:** パラメータ多様体 $\mathcal{M}$ 上の自然勾配法
3. **Fisher情報量幾何:** パラメータ空間の計量として Fisher情報行列を使用

$$g_{ij} = \sum_n \frac{\partial \mathcal{G}_0^{\text{And}}}{\partial \theta_i}^* \frac{\partial \mathcal{G}_0^{\text{And}}}{\partial \theta_j}$$

ここで $\theta = \{\epsilon_k, V_k\}$。この計量を用いた自然勾配降下は、パラメータの再パラメトライゼーションに対して不変であり、VA10AのHessian近似よりも幾何学的に自然な最適化を実現する。

### 4.3 Galerkin法に基づくモダン設計

Galerkin法の観点からのモジュール設計：

```
module green_function_space
  ! Green関数のHilbert空間を定義
  ! 内積、ノルム、射影演算子
end module

module anderson_model
  ! Anderson模型のパラメータ多様体
  ! G0_and の計算と、パラメータに関する微分
end module

module variational_solver
  ! Galerkin射影問題のソルバー
  ! 残差直交条件の反復解法
  ! 最適化アルゴリズムの抽象インターフェース
end module

module optimizer
  ! 最適化アルゴリズムの具体実装
  ! - BFGS法（VA10Aの等価実装）
  ! - L-BFGS法
  ! - 信頼領域法
  ! - 自然勾配法
end module
```

## 5. 学術的推奨手法：本問題に最適なアプローチ

本問題 — Anderson不純物模型パラメータのフィッティング — の特性を考慮し、学術的に最も適切な手法を推奨順に提示する。

### 5.1 第一推奨：Levenberg-Marquardt法（非線形最小二乗法）

**推奨理由:** 本問題は本質的に**非線形最小二乗問題**である。

目的関数は残差の和：

$$\mathcal{F}(\theta) = \sum_{n=0}^{N_\omega} |r_n(\theta)|^2, \quad r_n(\theta) = \mathcal{G}_0(i\omega_n) - \mathcal{G}_0^{\text{And}}(i\omega_n; \theta)$$

VA10Aのような一般的な最適化法はこの**和の構造を無視**している。Levenberg-Marquardt法はJacobian行列 $J_{nj} = \partial r_n / \partial \theta_j$ を用いてHessian近似 $B \approx J^T J$ を構成する。これは以下の利点をもたらす：

- Hessianの二次項 $\sum_n r_n \nabla^2 r_n$ を無視できる（残差が小さい場合に正当化される）
- $J^T J$ は常に半正定値 → Hessianの正定値性が自動的に保証される
- 信頼領域法と自然に統合される

**学術的根拠:**
- K. Levenberg, "A Method for the Solution of Certain Non-Linear Problems in Least Squares," Quarterly of Applied Mathematics 2, 164-168 (1944)
- D.W. Marquardt, "An Algorithm for Least-Squares Estimation of Nonlinear Parameters," SIAM J. Appl. Math. 11, 431-441 (1963)
- J.J. Moré, "The Levenberg-Marquardt Algorithm: Implementation and Theory," in Numerical Analysis, Lecture Notes in Mathematics 630, Springer (1978)

**DMFT文脈での実績:**
Levenberg-Marquardt法はDMFTのバスフィッティングにおいて広く使用されている。特にCT-QMC（連続時間量子モンテカルロ）法との組み合わせで標準的な手法となっている。

**実装:** MINPACK ライブラリの `lmder` / `lmdif` が標準的実装。Fortranで書かれており、直接利用可能。
- J.J. Moré, B.S. Garbow, K.E. Hillstrom, "User Guide for MINPACK-1," Argonne National Laboratory Report ANL-80-74 (1980)

### 5.2 第二推奨：Broyden法による自己無撞着方程式の直接求解

**推奨理由:** DMFTの自己無撞着条件を**不動点問題**として直接扱う。

最適化問題としてではなく、非線形方程式系：

$$\mathcal{G}_0^{\text{And}}(i\omega_n; \theta^*) = \mathcal{G}_0(i\omega_n) \quad \forall n$$

の根を求める問題として定式化する。Broydenの良い方法（Broyden's first method）は：

$$\theta_{k+1} = \theta_k - B_k^{-1} F(\theta_k)$$

ここで $F(\theta) = (\mathcal{G}_0^{\text{And}} - \mathcal{G}_0)$ であり、$B_k$ はJacobianの近似を秘割線条件で更新する。

**学術的根拠:**
- C.G. Broyden, "A Class of Methods for Solving Nonlinear Simultaneous Equations," Mathematics of Computation 19, 577-593 (1965)
- D.D. Johnson, "Modified Broyden's method for accelerating convergence in self-consistent calculations," Phys. Rev. B 38, 12807 (1988) — **DMFT分野で標準的に引用される**

**利点:**
- DMFTの自己無撞着ループの収束加速に直接適用可能
- Anderson mixing（線形混合）の一般化として理解できる
- $N_\omega \gg n_\text{param}$ の場合（本問題に該当）に効率的

### 5.3 第三推奨：Potthoffの自己エネルギー汎関数理論（SFT）

**推奨理由:** 最適化問題そのものを物理学の変分原理から再導出する。

Potthoff (2003) の自己エネルギー汎関数理論は、Luttinger-Ward汎関数の Legendre変換として大ポテンシャル：

$$\Omega[\Sigma] = F[\Sigma] + \text{Tr}\ln\left[-(\mathcal{G}_0^{-1} - \Sigma)^{-1}\right]$$

の停留条件 $\delta\Omega / \delta\Sigma = 0$ からDMFTの自己無撞着条件を導く。Anderson模型パラメータ空間 $\{\epsilon_k, V_k\}$ への制限は、この汎関数の**離散変分**を与える：

$$\frac{\partial \Omega}{\partial \epsilon_k} = 0, \quad \frac{\partial \Omega}{\partial V_k} = 0$$

**学術的根拠:**
- M. Potthoff, "Self-energy-functional approach to systems of correlated electrons," Eur. Phys. J. B 32, 429-436 (2003)
- M. Potthoff, "Self-energy-functional theory," in Strongly Correlated Systems, Springer Series in Solid-State Sciences 171 (2012)

**意義:** バスフィッティングの目的関数の選び方（L1ノルム vs L2ノルム、周波数の重み付け等）の恣意性を排除し、**熱力学的に整合した変分原理**から最適化問題を一意に導出できる。

### 5.4 第四推奨：有理関数近似（Padé近似 / AAA アルゴリズム）

**推奨理由:** Anderson模型Green関数のフィッティングは、本質的に**有理関数近似問題**である。

$\mathcal{G}_0^{\text{And}}(z)$ は $z$ の有理関数（部分分数展開）であるから、ターゲット関数 $\mathcal{G}_0(z)$ の最良有理近似を求める問題に帰着する：

$$\mathcal{G}_0(z) \approx \frac{p(z)}{q(z)} = \frac{a_0 + a_1 z + \cdots}{b_0 + b_1 z + \cdots + z^{N_s}}$$

**最新の学術的手法:**

1. **AAAアルゴリズム** (Adaptive Antoulas-Anderson):
   - Y. Nakatsukasa, O. Sète, L.N. Trefethen, "The AAA Algorithm for Rational Approximation," SIAM J. Sci. Comput. 40, A1494-A1522 (2018)
   - 貪欲アルゴリズムによるbarycentric有理近似。ロバストで最適化不要

2. **Vector Fitting (VF)**:
   - B. Gustavsen, A. Semlyen, "Rational Approximation of Frequency Domain Responses by Vector Fitting," IEEE Trans. Power Delivery 14, 1052-1061 (1999)
   - 電気工学で標準的。極と留数を反復的に求める

3. **Nevanlinna解析接続**:
   - J. Fei, C.-N. Yeh, E. Gull, "Nevanlinna Analytical Continuation," Phys. Rev. Lett. 126, 056402 (2021)
   - 因果律（正のスペクトル関数）を厳密に保証する有理関数近似

**意義:** 最適化問題を回避し、有理関数近似の数学的理論（Padé近似の一意性、最良有理近似のde la Vallée-Poussinの定理）に基づくアプローチ。特にAAAアルゴリズムは反復最適化なしに高精度な有理近似を与える。

### 5.5 推奨手法の比較

| 手法 | 理論的基盤 | 実装難易度 | VA10Aとの互換性 | 物理的整合性 |
|------|-----------|-----------|----------------|-------------|
| Levenberg-Marquardt | 非線形最小二乗 | 低（MINPACK利用） | 高（目的関数同一） | 中 |
| Broyden法 | 非線形方程式 | 中 | 中（定式化変更） | 高 |
| SFT変分法 | 熱力学的変分原理 | 高 | 低（目的関数変更） | 最高 |
| AAA/有理近似 | 近似理論 | 中 | 低（アプローチ変更） | 高 |

### 5.6 実践的推奨

**段階1（即座に実装可能）:** Levenberg-Marquardt法（MINPACK）への置き換え。目的関数の構造を活用でき、VA10Aと同じ入出力で動作させられる。残差ベクトル $r_n$ を個別に返すようにenergyサブルーチンを修正するだけでよい。

**段階2（中期）:** Broyden法によるDMFT自己無撞着ループ全体の加速。minimizeの置き換えに留まらず、外側の反復ループの収束性を改善する。

**段階3（長期・研究指向）:** SFT変分法またはAAA有理近似による問題の再定式化。物理的により深い理解に基づくアプローチで、バスフィッティングの精度と一意性を根本的に改善する。

## 6. 等価プログラム構築への道筋

### 6.1 段階的移行戦略

1. **第1段階（等価性保証）:** VA10Aを現代Fortranで忠実に再実装し、テストで等価性を検証
2. **第2段階（モジュール化）:** 最適化器、目的関数、Green関数計算を独立モジュールに分離
3. **第3段階（理論的拡張）:** Galerkin法の枠組みで抽象インターフェースを設計し、複数の最適化アルゴリズムを差し替え可能にする
4. **第4段階（物理的拡張）:** 自然勾配法やFisher情報計量の導入で、物理的に意味のある最適化を実現

### 6.2 現代Fortranの活用

| 旧来の手法 | 現代的手法 |
|-----------|-----------|
| COMMON block | MODULE変数 |
| EXTERNAL宣言 | ABSTRACT INTERFACE |
| GOTO文による制御フロー | DO/SELECT CASE/EXIT/CYCLE |
| パック行列格納 | 派生型(TYPE)によるカプセル化 |
| implicit型宣言 | IMPLICIT NONE + 明示的型宣言 |
| 作業配列w(*) | ALLOCATABLE配列 |
| 数値微分 | 自動微分またはユーザ提供勾配 |

### 6.3 検証基準

モダナイゼーション後のコードは、以下を満たすべきである：

1. **数値的等価性:** 同一の入力に対して、オリジナルと同一の最適化結果（浮動小数点精度内）
2. **テスト互換性:** Rosenbrock関数と二次関数のテストが同一のPASS結果
3. **理論的正当性:** Galerkin射影の残差直交条件が満たされていることの直接検証
4. **拡張性:** 新しい最適化アルゴリズムの追加がインターフェース変更なしに可能

## 7. 結論

lisalanc.fの`minimize`（VA10A）は、準Newton共役勾配法のBroyden族更新をCholesky分解上で実装した高度なアルゴリズムである。その物理学的文脈において、これはGreen関数空間上のGalerkin射影問題を解くものとして再解釈できる。

真のモダナイゼーションとは、この理論的構造を明示化し、現代的な数値最適化理論（L-BFGS、信頼領域法、自動微分、自然勾配法）と関数解析的枠組み（Hilbert空間上の変分問題、Riemannian最適化）を統合した設計を実現することである。
