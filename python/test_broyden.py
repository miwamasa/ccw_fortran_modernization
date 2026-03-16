"""
Broyden法 Python版のテスト

Fortran版と同一のテスト関数で検証し、結果を比較する。

Test 1: Rosenbrock関数 (残差形式)
    r1 = 1 - x,  r2 = 10*(y - x^2)
    f = r1^2 + r2^2 = (1-x)^2 + 100*(y-x^2)^2
    最小値: f(1, 1) = 0

Test 2: 二次関数 (残差形式)
    r1 = x - 3,  r2 = y + 2
    f = r1^2 + r2^2 = (x-3)^2 + (y+2)^2
    最小値: f(3, -2) = 0

Test 3: 3変数連立方程式（線形系）
    r1 = 2x - y + z - 5
    r2 = x + 3y - z - 2
    r3 = -x + y + 2z - 11
    解: (1, 2, 5)
"""

import numpy as np
import sys
import os

# broyden_solver.py を同ディレクトリからインポート
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from broyden_solver import broyden_minimize, broyden_solve


def rosenbrock_resid(x):
    """Rosenbrock関数の残差: r1 = 1-x, r2 = 10*(y-x^2)"""
    return np.array([1.0 - x[0], 10.0 * (x[1] - x[0]**2)])


def quadratic_resid(x):
    """二次関数の残差: r1 = x-3, r2 = y+2"""
    return np.array([x[0] - 3.0, x[1] + 2.0])


def system3_resid(x):
    """3変数連立方程式（一意解を持つ線形寄りの系）
    r1 = 2*x - y + z - 5 = 0  → x=1: 2-2+1-5=-4 NG
    解: x=1, y=2, z=5  → 2-2+5-5=0 OK
    r2 = x + 3*y - z - 2 = 0  → 1+6-5-2=0 OK
    r3 = -x + y + 2*z - 11 = 0 → -1+2+10-11=0 OK
    """
    return np.array([
        2.0*x[0] - x[1] + x[2] - 5.0,
        x[0] + 3.0*x[1] - x[2] - 2.0,
        -x[0] + x[1] + 2.0*x[2] - 11.0,
    ])


def run_tests():
    npass = 0
    nfail = 0

    print("=" * 60)
    print("Broyden solver (Python) test suite")
    print("=" * 60)
    print()

    # ------------------------------------------------------------------
    # Test 1: Rosenbrock
    # ------------------------------------------------------------------
    print("Test 1: Rosenbrock function (Python Broyden)")
    print("  r1 = 1-x, r2 = 10*(y-x^2)")
    print("  Expected solution: (1, 1), f = 0")
    print()

    x0 = np.array([-1.0, 1.0])
    result = broyden_minimize(rosenbrock_resid, x0, eps=1e-12, maxiter=10000)

    print(f"  Result: x = {result['x'][0]:12.6f}, y = {result['x'][1]:12.6f}")
    print(f"  f(x,y) = {result['f']:.7e}")
    print(f"  iexit  = {result['iexit']}")
    print(f"  niter  = {result['niter']}")
    print(f"  nfev   = {result['nfev']}")
    print()

    tol = 1e-4
    err_x = abs(result['x'][0] - 1.0)
    err_y = abs(result['x'][1] - 1.0)
    if err_x < tol and err_y < tol:
        print("  >>> PASS: Converged to correct solution")
        npass += 1
    else:
        print(f"  >>> FAIL: |x-1| = {err_x:.4e}, |y-1| = {err_y:.4e}")
        nfail += 1
    print()

    # ------------------------------------------------------------------
    # Test 2: Quadratic
    # ------------------------------------------------------------------
    print("Test 2: Quadratic function (Python Broyden)")
    print("  r1 = x-3, r2 = y+2")
    print("  Expected solution: (3, -2), f = 0")
    print()

    x0 = np.array([0.0, 0.0])
    result = broyden_minimize(quadratic_resid, x0, eps=1e-12, maxiter=10000)

    print(f"  Result: x = {result['x'][0]:12.6f}, y = {result['x'][1]:12.6f}")
    print(f"  f(x,y) = {result['f']:.7e}")
    print(f"  iexit  = {result['iexit']}")
    print(f"  niter  = {result['niter']}")
    print(f"  nfev   = {result['nfev']}")
    print()

    tol = 1e-6
    err_x = abs(result['x'][0] - 3.0)
    err_y = abs(result['x'][1] + 2.0)
    if err_x < tol and err_y < tol:
        print("  >>> PASS: Converged to correct solution")
        npass += 1
    else:
        print(f"  >>> FAIL: |x-3| = {err_x:.4e}, |y+2| = {err_y:.4e}")
        nfail += 1
    print()

    # ------------------------------------------------------------------
    # Test 3: 3-variable system
    # ------------------------------------------------------------------
    print("Test 3: 3-variable linear system (Python Broyden)")
    print("  r1 = 2x-y+z-5, r2 = x+3y-z-2, r3 = -x+y+2z-11")
    print("  Expected solution: (1, 2, 5)")
    print()

    x0 = np.array([0.0, 0.0, 0.0])
    result = broyden_solve(system3_resid, x0, eps=1e-10, maxiter=10000)

    print(f"  Result: x = {result['x'][0]:12.6f}, y = {result['x'][1]:12.6f}, z = {result['x'][2]:12.6f}")
    print(f"  ||F||  = {result['fnorm']:.7e}")
    print(f"  iexit  = {result['iexit']}")
    print(f"  niter  = {result['niter']}")
    print(f"  nfev   = {result['nfev']}")
    print()

    tol = 1e-6
    expected = np.array([1.0, 2.0, 5.0])
    errors = np.abs(result['x'] - expected)
    if np.all(errors < tol):
        print("  >>> PASS: Converged to correct solution")
        npass += 1
    else:
        print(f"  >>> FAIL: errors = {errors}")
        nfail += 1
    print()

    # ------------------------------------------------------------------
    # Fortran版との比較
    # ------------------------------------------------------------------
    print("=" * 60)
    print("Fortran版との比較（Rosenbrock関数）")
    print("=" * 60)
    print()
    print("  | 実装      | f(x,y)       | iexit |")
    print("  |-----------|-------------|-------|")
    print(f"  | Fortran   | 4.93e-30    | 1     |")

    x0 = np.array([-1.0, 1.0])
    result = broyden_minimize(rosenbrock_resid, x0, eps=1e-12)
    print(f"  | Python    | {result['f']:.2e}    | {result['iexit']}     |")
    print()

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    print("=" * 60)
    print(f"  Results: {npass} PASS, {nfail} FAIL")
    print("=" * 60)

    return nfail


if __name__ == "__main__":
    nfail = run_tests()
    sys.exit(1 if nfail > 0 else 0)
