# lisalanc.f 分析レポート

## 1. プログラム全体の概要

`lisalanc.f`は **LISA (Local Impurity Self Consistent Approximation)** プログラムであり、強相関フェルミオン系（特に単一バンドHubbardモデル）の自己無撞着計算を行う研究コードである。

- **著者:** W. Krauth, M. Caffarel
- **発表:** A. Georges, G. Kotliar, W. Krauth, M. Rozenberg, Reviews of Modern Physics (1996)
- **アルゴリズム:** M. Caffarel and W. Krauth, Phys. Rev. Lett. 72, 1545 (1994)

### 物理的背景

Hubbardモデルの常磁性相における零温度の熱力学的Green関数を計算する。Anderson不純物モデルのLanczos対角化による解法を用い、動的平均場理論（DMFT）の自己無撞着ループを実行する。

### 全体フロー

```
main (lisalanc)
  ├── initial: パラメータ初期化、andparファイル読み込み
  ├── diag: Lanczos対角化
  │    ├── setstate: Hilbert空間の構築
  │    ├── findgroundstatel: Lanczos法＋QL対角化で基底状態計算
  │    │    ├── lanczos: Lanczos反復
  │    │    ├── tql2: QL法による固有値計算（SLATEC）
  │    │    └── findgroundstate: 射影法による精密化
  │    ├── computegreen: Green関数計算
  │    └── amultv: 疎行列-ベクトル積
  ├── search: Anderson模型パラメータ最適化
  │    ├── minimize: 準Newton共役勾配法（VA10A）
  │    └── energy: 目的関数（G0の不一致度）
  └── calcg0: Anderson模型Green関数計算
```

## 2. 主要サブルーチン

| サブルーチン | 行番号 | 目的 |
|-------------|--------|------|
| `lisalanc` (main) | 61-192 | メインプログラム。入力読み込み、対角化、自己無撞着ループ |
| `initial` | 200-215 | Anderson模型パラメータの初期化 |
| `rheader`/`wheader` | 217-256 | andparファイルの読み書き |
| `diag` | 271-420 | Lanczos対角化。各粒子数セクターで基底状態を計算 |
| `setstate` | 430-560 | Hilbert空間の構築、Hamiltonian行列要素の計算 |
| `search` | 707-746 | minimize を呼び出してAnderson模型パラメータを最適化 |
| `computegreen` | 758-888 | Lanczos法によるGreen関数への寄与を計算 |
| `minimize` | 890-1193 | 準Newton共役勾配法による最適化（詳細は下記） |
| `energy` | 1203-1224 | 目的関数: G0_andとG0_newの不一致度 |
| `calcg0` | 1233-1248 | Anderson模型のGreen関数G0を計算 |
| `amultv` | 1258-1266 | 疎Hamiltonian行列とベクトルの積 |
| `findgroundstatel` | 1268-1370 | Lanczos＋QL法で基底状態を求める |
| `lanczos` | 1375-1420 | 単一Lanczos反復ステップ |
| `findgroundstate` | 1425-1456 | 射影法による基底状態の精密化 |
| `tql2` | 1460-1620 | QL法による三重対角行列の固有値計算（SLATEC） |

## 3. minimizeサブルーチンの詳細分析

### 3.1 概要

```fortran
subroutine minimize(funct, n, x, f, g, h, w, dfn, xm, hh, eps, mode, maxfn, iprint, iexit)
```

- **行番号:** 890-1193
- **アルゴリズム:** 準Newton共役勾配法
- **原名:** VA10A（Harwellライブラリ由来と推定）
- **特徴:** コメントに「最も信頼できる共役勾配ルーチン」と記載。原作者不明。

### 3.2 引数一覧

| # | 引数 | 型 | 入出力 | 説明 |
|---|------|-----|--------|------|
| 1 | `funct` | EXTERNAL | 入力 | 目的関数。`call funct(n, x, f)` の形式で呼ばれる |
| 2 | `n` | INTEGER | 入力 | 最適化変数の数 |
| 3 | `x(n)` | DOUBLE PRECISION | 入出力 | 入力：初期推定値、出力：最適化された解 |
| 4 | `f` | DOUBLE PRECISION | 出力 | 最小化された関数値 |
| 5 | `g(n)` | DOUBLE PRECISION | 作業用 | 勾配ベクトル（数値微分で計算） |
| 6 | `h(n*(n+1)/2)` | DOUBLE PRECISION | 入出力 | Hessian行列（対称行列のパック形式、下三角） |
| 7 | `w(*)` | DOUBLE PRECISION | 作業用 | 作業配列（サイズ ≥ n*(n+4)程度） |
| 8 | `dfn` | DOUBLE PRECISION | 入力 | 期待される関数値の減少量。ステップサイズの初期推定に使用 |
| 9 | `xm(n)` | DOUBLE PRECISION | 入力 | 各変数のスケール因子／最大ステップ幅 |
| 10 | `hh` | DOUBLE PRECISION | 入力 | 数値微分のステップ幅（通常1.0e-5） |
| 11 | `eps` | DOUBLE PRECISION | 入力 | 収束判定の許容誤差（通常0.00001） |
| 12 | `mode` | INTEGER | 入力 | 初期化モード: 1=単位行列、2=与えられたHessian使用、3=ユーザ指定状態 |
| 13 | `maxfn` | INTEGER | 入力 | 関数評価の最大回数 |
| 14 | `iprint` | INTEGER | 入力 | 出力制御: 0=出力なし、>0=iprint反復ごとに出力 |
| 15 | `iexit` | INTEGER | 出力 | 終了状態: 0=初期、1=収束、2=勾配条件不満足、3=最大回数到達 |

### 3.3 アルゴリズムの詳細

#### フェーズ1: 初期化（行919-963）

- **mode=1:** Hessian行列を単位行列で初期化
- **mode=2:** 与えられたHessian行列のCholesky分解を検証
- **mode=3:** ユーザが提供する状態をそのまま使用
- Hessianの最小対角要素`dmin`を計算、正定値性を確認

#### フェーズ2: 勾配計算（行1171-1192）

2つのモード:
- **前進差分** (idiff=1): `g(i) = (f(x+h) - f(x)) / h`
- **中心差分** (idiff=2): `g(i) = (f(x+h) - f(x-h)) / (2h)`

ステップ幅: `z = hh * xm(i)`（各変数のスケールに応じた幅）

最初は前進差分で計算し、収束後に中心差分に切り替えて精度を上げる（行1160-1162）。

#### フェーズ3: 探索方向の計算（行988-1016）

Hessianの逆行列と勾配から探索方向を計算:
1. `w = -H⁻¹ * g`（Cholesky分解を用いた三角行列の前進・後退代入）
2. 探索方向のスケーリング確認: `z = max(|w(i)|/xm(i))`
3. 勾配と探索方向の内積: `gs0 = g·w`（下降方向であることを確認）
4. 初期ステップ幅: `alpha = min(1, -2*df/gs0)`

#### フェーズ4: 線探索（行1026-1082）

適応的な線探索:
1. **ステップ倍増**（行1036-1054）: f1 < fの場合、alphaを2倍にして探索続行
2. **ステップ半減**（行1055-1078）: f1 ≥ fの場合、alphaを半分にして補間
3. **二次補間**（行1072-1074）: f, f1, f2から放物線近似で最適ステップ推定
4. 収束判定: `tot < aeps`（ステップが十分小さい場合）

#### フェーズ5: Hessian更新（行1096-1154）

Broyden族の更新公式:
1. 勾配の変化量から更新ベクトルを計算
2. Cholesky因子を直接更新（行列の再分解を回避）
3. Hessianの正定値性を維持（`z ≤ 0`の場合は`dmin`で置換）

2回の更新サイクル（link=1, link=2）でBroyden族のランク2更新を実現。

#### フェーズ6: 収束判定

- **iexit=1:** 収束。`tot < aeps`（ステップ幅が許容誤差以下）
- **iexit=2:** 勾配条件不満足。`gs0 ≥ 0`（探索方向が下降方向でない）
- **iexit=3:** 最大関数評価回数到達。`ifn ≥ maxfn`

### 3.4 searchサブルーチンからの呼び出し

```fortran
subroutine search(fmin, Nitermax)
  ! パラメータ数: nbparm = 2*Ns - 2
  !   Ns-1個のエネルギー準位 Epsk(2),...,Epsk(Ns)
  !   Ns-1個の結合強度 V(1),...,V(Ns-1)

  ! スケール因子: 各パラメータの絶対値 + 1e-15
  xprmt(i) = dabs(xtemp(i)) + 1.d-15

  ! 呼び出し
  call minimize(energy, nbparm, xtemp, fmin, g, hess, w,
                dfn=-0.5, xprmt, hh=1.e-5, deps=0.00001,
                mode=1, Nitermax=4000, iprint=0, iexit)
```

目的関数`energy`はAnderson模型のGreen関数G0_andと、Hubbardモデルから計算したG0_newの松原周波数上での不一致度を計算する:

```fortran
f = (1/(Iwmax+1)) * Σ |G0w(i) - G0wand(i)|
```

### 3.5 コードの特徴と問題点

- **GOTO文の多用:** 約30箇所のGOTO文がフロー制御に使用されている
- **計算GOTO文:** `go to (18, 54), link` や `go to (60, 20), link` のような計算GOTOが状態遷移に使用
- **暗黙の型宣言:** `implicit double precision (a-h,o-z)`による暗黙的型付け
- **パック行列格納:** Hessian行列が`h(n*(n+1)/2)`のパック形式で格納
- **マジックナンバー:** `7.0d0`, `5.0d0`, `12.0d0`, `0.1d0`等の係数が説明なしに使用

## 4. データ構造（lisalanc.dat）

```fortran
parameter(Nmpara=100)           ! 最適化パラメータの最大数
parameter(Iwmax=2**6)           ! 松原周波数の最大点数（64）
parameter(Nl=100, Tiny=1.e-13)  ! Lancozベクトル数、微小量
parameter(Ns=4, Nslmax=70, Ncontx=17) ! サイト数、Hilbert空間サイズ
parameter(Nitmax=1000, Ip=2*Ns) ! 最大反復数、全粒子数

! Green関数
complex*16 Gw(0:Iwmax)          ! 熱力学的Green関数
complex*16 G0w(0:Iwmax)         ! 裸のGreen関数
complex*16 G0wand(0:Iwmax)      ! Anderson模型のGreen関数

! Anderson模型パラメータ
dimension V(Ns)                  ! 結合強度
dimension Epsk(Ns)               ! エネルギー準位
```

## 5. 入出力ファイル

| ファイル | 種類 | 内容 |
|---------|------|------|
| `lisalanc.input` | 入力 | Beta（逆温度）、U（相互作用）、xmu（化学ポテンシャル）、Cutoff |
| `lisalanc.andpar` | 入出力 | Anderson模型パラメータ（Epsk, V の初期値と更新値） |
| `lisalanc.green` | 出力 | Green関数 G(iω_n) |
| `lisalanc.diff` | 出力 | G0_wとG0_wandの比較 |
