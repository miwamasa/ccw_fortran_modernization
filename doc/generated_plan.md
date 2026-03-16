# Fortranモダナイゼーション実装計画

## Context

`lisalanc.f`はLISA (Local Impurity Self Consistent Approximation) プログラムで、強相関フェルミオン系の自己無撞着計算を行う。特に`minimize`サブルーチン（890-1193行目）は、VA10Aとして知られる準Newton共役勾配法の最適化ルーチンで、GOTO文を多用した非常に難解なコードである。instructions.mdの4つのタスクを順番に実行する。

## タスク一覧

### タスク1: 分析 → `doc/analysis.md`

lisalanc.fの全体構造とminimizeルーチンの詳細分析を文書化する。

**内容:**
- プログラム全体の目的：Hubbardモデルの自己無撞着Lanczos対角化
- 主要サブルーチン構成（main, diag, search, computegreen, minimize, energy, calcg0, amultv, findgroundstatel, lanczos, tql2）
- **minimize詳細分析:**
  - アルゴリズム：準Newton共役勾配法（VA10A）
  - 数値勾配（前進差分/中心差分）、Broyden型Hessian更新、適応的線探索
  - 全15引数の意味と型
  - 収束判定・終了条件（iexit=1,2,3）
  - searchサブルーチンからの呼び出しパターン

### タスク2: テスト生成 → `test/`

minimizeを単体で取り出し、解析解が既知の2つの非線形関数でテストする。

**ファイル構成:**
- `test/minimize.f` — minimizeサブルーチンを抽出（lisalanc.fの890-1193行をそのままコピー）
- `test/test_minimize.f` — テストドライバープログラム
- `test/CMakeLists.txt` — CMakeビルドファイル
- ルート `CMakeLists.txt` — トップレベルCMakeファイル

**テスト関数（2つ）:**

1. **Rosenbrock関数** (2変数): f(x,y) = (1-x)² + 100(y-x²)²
   - 最小値: f(1,1) = 0
   - 非線形で共役勾配法のよいテスト

2. **二次関数** (2変数): f(x,y) = (x-3)² + (y+2)²
   - 最小値: f(3,-2) = 0
   - 簡単な凸関数で基本動作確認

**テストドライバーの実装:**
```
program test_minimize
  - Rosenbrock関数テスト: 初期点(-1,1)から開始、(1,1)への収束を確認
  - 二次関数テスト: 初期点(0,0)から開始、(3,-2)への収束を確認
  - 許容誤差内（1e-3程度）で解に到達しているか判定
  - PASS/FAILを出力
end program
```

**CMakeLists.txt:**
```cmake
cmake_minimum_required(VERSION 3.10)
project(minimize_test Fortran)
enable_testing()
add_executable(test_minimize test/test_minimize.f test/minimize.f)
add_test(NAME minimize_test COMMAND test_minimize)
```

### タスク3: 初期テスト実行 → `doc/original_test_result.md`

- CMakeでビルド・実行
- テスト結果（各関数の収束結果、iexit値、関数評価回数）を記録
- `doc/original_test_result.md`に結果を文書化

### タスク4: 理論構築 → `doc/theory.md`

**目的:** 単なるGOTO排除リファクタリングではなく、高い視座からの理論的背景を構築

**構成:**
1. **元のアルゴリズムの理論的位置づけ**
   - VA10Aの数学的基盤：準Newton法、Broyden族の更新式
   - 共役勾配法の理論（Fletcher-Reeves, Polak-Ribière）

2. **物理学的文脈での最適化理論**
   - Anderson不純物模型のパラメータフィッティングとしての定式化
   - Green関数の汎関数としての変分原理
   - DMFT (動的平均場理論) における自己無撞着条件の変分的解釈

3. **モダナイゼーションの理論的基盤**
   - **変分法的アプローチ:** 目的関数をGreen関数空間上の汎関数として定式化
   - **Galerkin法との対応:** Anderson模型パラメータ空間への射影としてのGalerkin的解釈
   - **現代的最適化フレームワーク:** L-BFGS法、信頼領域法、自動微分との関係
   - **関数空間上の勾配降下:** Hilbert空間での最急降下法としての再定式化

4. **等価プログラム構築への道筋**
   - モジュール化設計（最適化器、目的関数、Green関数計算の分離）
   - 現代Fortranの機能活用（モジュール、派生型、インターフェース）
   - 抽象的な最適化フレームワークの設計

## 修正対象ファイル

| ファイル | 操作 |
|---------|------|
| `doc/analysis.md` | 新規作成 |
| `test/minimize.f` | 新規作成（lisalanc.fから抽出） |
| `test/test_minimize.f` | 新規作成 |
| `test/CMakeLists.txt` | 新規作成 |
| `CMakeLists.txt` | 新規作成（トップレベル） |
| `doc/original_test_result.md` | 新規作成（テスト結果） |
| `doc/theory.md` | 新規作成 |

## 検証方法

1. `mkdir build && cd build && cmake .. && make` でビルドが通ること
2. `ctest` または `./test_minimize` でテストが実行され、両関数がPASSすること
3. doc/以下の3つの文書が適切な内容で生成されていること
