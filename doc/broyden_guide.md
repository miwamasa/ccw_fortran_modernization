# Broyden法の理論と実装ガイド

## 1. 概要

Broyden法は、非線形方程式系 $F(x) = 0$ を解くための準Newton法である。Newton法がステップごとにJacobian行列を再計算するのに対し、Broyden法は初期Jacobianからの**ランク1更新**によって近似Jacobianを効率的に維持する。

- **提案者:** C.G. Broyden (1965)
- **論文:** "A Class of Methods for Solving Nonlinear Simultaneous Equations," Mathematics of Computation 19, 577-593
- **DMFT応用:** D.D. Johnson, Phys. Rev. B 38, 12807 (1988)

## 2. 数学的基礎

### 2.1 Newton法からの出発

$n$ 変数の非線形方程式系 $F: \mathbb{R}^n \to \mathbb{R}^n$ に対して、Newton法は：

$$x_{k+1} = x_k - J(x_k)^{-1} F(x_k)$$

ここで $J(x_k) = \frac{\partial F}{\partial x}\bigg|_{x_k}$ はJacobian行列。

**Newton法の問題点:**
- 毎ステップで $n^2$ 個の偏微分の計算が必要（$n+1$ 回の関数評価）
- $n$ が大きい場合に計算コストが高い

### 2.2 Broyden法の着想

**秘割線条件 (Secant condition):** 近似Jacobian $B_k$ が以下を満たすことを要求する：

$$B_{k+1} s_k = y_k$$

ここで：
- $s_k = x_{k+1} - x_k$（ステップ）
- $y_k = F(x_{k+1}) - F(x_k)$（関数値の変化）

これは $n$ 本の方程式に $n^2$ 個の未知数（$B_{k+1}$ の成分）があるため、不定系である。Broydenは**最小変更条件**を追加する：

$$B_{k+1} = \arg\min_B \| B - B_k \|_F \quad \text{subject to} \quad B s_k = y_k$$

### 2.3 Broyden更新公式

上記の条件付き最適化を解くと、**Broyden's first method**（"good" Broyden）が得られる：

$$B_{k+1} = B_k + \frac{(y_k - B_k s_k) s_k^T}{s_k^T s_k}$$

これはランク1更新であり、$O(n^2)$ の計算量でJacobian近似を更新できる。

### 2.4 逆行列の直接更新（Sherman-Morrison公式）

実用上は $B_k^{-1}$ を直接保持し、**Sherman-Morrison公式**で更新する：

$$B_{k+1}^{-1} = B_k^{-1} + \frac{(s_k - B_k^{-1} y_k) s_k^T B_k^{-1}}{s_k^T B_k^{-1} y_k}$$

**利点:**
- 連立方程式を解く必要がない（$B^{-1}F$ を直接計算）
- 更新も $O(n^2)$ の外積演算のみ
- 数値的に安定

### 2.5 アルゴリズム

```
Input: F(x), x₀, ε (許容誤差), maxiter
Output: x* (F(x*)≈0 の解)

1. B₀⁻¹ ← J(x₀)⁻¹  (初期Jacobianの逆行列を数値微分で計算)
2. F₀ ← F(x₀)
3. for k = 0, 1, 2, ... do
4.   dx ← -B_k⁻¹ F_k                     [Newton-like ステップ]
5.   x_{k+1} ← x_k + dx
6.   F_{k+1} ← F(x_{k+1})
7.   if ||F_{k+1}|| < ε then return x_{k+1}  [収束判定]
8.   s ← dx
9.   y ← F_{k+1} - F_k
10.  B_{k+1}⁻¹ ← B_k⁻¹ + (s - B_k⁻¹ y)(s^T B_k⁻¹) / (s^T B_k⁻¹ y)
                                            [Sherman-Morrison更新]
11. end for
```

## 3. 収束性

### 3.1 局所的超線形収束

初期点が解に十分近く、初期Jacobian近似 $B_0$ が十分良い場合、Broyden法は**超線形収束 (superlinear convergence)** を示す：

$$\lim_{k \to \infty} \frac{\|x_{k+1} - x^{*}\|}{\|x_k - x^{*}\|} = 0$$

これはNewton法の二次収束より遅いが、線形収束より速い。

**根拠:** Dennis-Moré条件 (1974) — 秘割線条件を満たすBroyden更新は、Dennis-Moré条件を満たし、超線形収束を保証する。

### 3.2 大域的収束

Broyden法は一般に大域的収束を保証しない。実用上は以下の手法で安定化する：

1. **ステップ制限:** $\|dx\| > \Delta_{\max}$ の場合にスケーリング
2. **線探索:** $\|F(x_k + \alpha dx)\|$ を最小化する $\alpha$ を探索
3. **再初期化:** 収束が停滞した場合にJacobianを再計算

## 4. DMFT文脈での応用

### 4.1 バスフィッティングへの適用

DMFTの自己無撞着条件は、Anderson模型パラメータ $\theta = \{\epsilon_k, V_k\}$ に対する非線形方程式系：

$$F_n(\theta) = \mathcal{G}_0(i\omega_n) - \mathcal{G}_0^{\text{And}}(i\omega_n; \theta) = 0 \quad (n = 0, 1, \ldots, N_\omega)$$

$N_\omega > n_{\text{param}}$ の場合（過剰決定系）、最小二乗の意味で解く。$N_\omega = n_{\text{param}}$ の場合は直接Broyden法が適用可能。

### 4.2 Johnson mixing との関係

D.D. Johnson (1988) のDMFT向け手法は、Broyden法の変種で以下の特徴を持つ：

- 自己エネルギー $\Sigma(i\omega_n)$ に対する不動点反復の加速
- 過去 $m$ ステップの情報を保持する**limited-memory**版
- Anderson mixing の一般化として理解可能

$$\Sigma_{\text{new}} = \Sigma_{\text{old}} + \alpha \Delta\Sigma + \text{Broyden correction}$$

### 4.3 VA10A (minimize) との対応

| 観点 | VA10A | Broyden法 |
|------|-------|-----------|
| 解く問題 | $\min f(x)$ | $F(x) = 0$ |
| 更新対象 | Hessian $\nabla^2 f$ の近似 | Jacobian $\partial F / \partial x$ の近似 |
| 更新公式 | Broyden族ランク2更新 | Broyden ランク1更新 |
| 線探索 | 必要（ステップ幅 $\alpha$ を決定） | 通常不要（全ステップ採用） |
| 関数評価 | 毎ステップ $n$ 回（勾配用）+ 線探索 | 毎ステップ1回のみ |

## 5. 実装上の注意点

### 5.1 初期Jacobianの計算

初期Jacobianは数値微分で計算する：

$$J_{ij} \approx \frac{F_i(x + h e_j) - F_i(x)}{h}$$

ステップ幅 $h$ は $h \approx \sqrt{\epsilon_{\text{mach}}} \approx 10^{-8}$ (倍精度) が適切。

### 5.2 特異性の回避

Sherman-Morrison更新の分母 $s^T B^{-1} y$ がゼロに近い場合、更新をスキップする：

```python
if abs(s @ Binv_y) < 1e-30:
    continue  # 更新をスキップ
```

### 5.3 ステップ制限

発散防止のため、ステップサイズに上限を設ける：

```python
if norm(dx) > max_step:
    dx *= max_step / norm(dx)
```

## 6. Broyden法の変種

| 変種 | 更新公式 | 特徴 |
|------|---------|------|
| Broyden's first (good) | $B + \frac{(y-Bs)s^T}{s^Ts}$ | 標準的。$B$ を更新 |
| Broyden's second (bad) | $B + \frac{(y-Bs)s^TB^{-1}}{s^TB^{-1}y}$ | $B^{-1}$ を直接更新（本実装） |
| Anderson mixing | 過去 $m$ ステップの線形結合 | DMFT標準。limited-memory |
| Modified Broyden | Johnson (1988) | DMFT自己無撞着計算に最適化 |

## 7. 参考文献

1. C.G. Broyden, "A Class of Methods for Solving Nonlinear Simultaneous Equations," Math. Comp. 19, 577 (1965)
2. J.E. Dennis Jr., J.J. Moré, "A Characterization of Superlinear Convergence and Its Application to Quasi-Newton Methods," Math. Comp. 28, 549 (1974)
3. D.D. Johnson, "Modified Broyden's method for accelerating convergence in self-consistent calculations," Phys. Rev. B 38, 12807 (1988)
4. J. Nocedal, S.J. Wright, "Numerical Optimization," 2nd ed., Springer (2006), Chapter 11
5. A. Georges et al., "Dynamical mean-field theory of strongly correlated fermion systems and the limit of infinite dimensions," Rev. Mod. Phys. 68, 13 (1996)
