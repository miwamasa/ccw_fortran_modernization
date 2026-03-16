"""
Broyden法による非線形方程式系ソルバー

Broyden's first method ("good" Broyden) の実装。
初期JacobianをSherman-Morrison公式でランク1更新する。

参考文献:
    C.G. Broyden, Math. Comp. 19, 577 (1965)
    D.D. Johnson, Phys. Rev. B 38, 12807 (1988)
"""

import numpy as np
from typing import Callable, Optional


def broyden_solve(
    func: Callable[[np.ndarray], np.ndarray],
    x0: np.ndarray,
    eps: float = 1e-10,
    maxiter: int = 10000,
    max_step: float = 10.0,
    verbose: bool = False,
) -> dict:
    """
    Broyden法で非線形方程式系 F(x) = 0 を解く。

    Parameters
    ----------
    func : callable
        F(x) を返す関数。x: ndarray(n,) -> F: ndarray(n,)
    x0 : ndarray
        初期推定値
    eps : float
        収束判定許容誤差 (||F(x)|| < eps で収束)
    maxiter : int
        最大反復回数
    max_step : float
        ステップサイズの上限
    verbose : bool
        反復情報を出力するか

    Returns
    -------
    dict
        'x': 解, 'fnorm': 最終残差ノルム, 'niter': 反復回数,
        'nfev': 関数評価回数, 'iexit': 終了状態 (1=収束, 2=ステップ微小, 3=最大回数)
    """
    x = np.array(x0, dtype=float)
    n = len(x)
    nfev = 0

    # 初期Jacobianを数値微分で計算
    hh = 1e-7
    F0 = func(x)
    nfev += 1
    Jac = np.zeros((n, n))
    for j in range(n):
        x_pert = x.copy()
        x_pert[j] += hh
        F1 = func(x_pert)
        nfev += 1
        Jac[:, j] = (F1 - F0) / hh

    # Jacobian逆行列
    try:
        Binv = np.linalg.inv(Jac)
    except np.linalg.LinAlgError:
        Binv = np.eye(n)

    F_old = F0
    fnorm = np.linalg.norm(F_old)

    if verbose:
        print(f"Broyden solver: iter=0, ||F||={fnorm:.7e}")

    iexit = 0
    for k in range(1, maxiter + 1):
        # Newton-like ステップ
        dx = -Binv @ F_old

        # ステップ制限
        step_norm = np.linalg.norm(dx)
        if step_norm > max_step:
            dx *= max_step / step_norm
            step_norm = max_step

        # 更新
        x_new = x + dx
        F_new = func(x_new)
        nfev += 1
        fnorm = np.linalg.norm(F_new)

        if verbose:
            print(f"  iter={k}, ||F||={fnorm:.7e}, ||dx||={step_norm:.7e}")

        # 収束判定
        if fnorm < eps:
            x = x_new
            iexit = 1
            break

        if step_norm < eps * (np.linalg.norm(x) + eps):
            x = x_new
            iexit = 2
            break

        # Broyden逆行列更新 (Sherman-Morrison)
        s = dx
        y = F_new - F_old
        Binv_y = Binv @ y
        sBinv = s @ Binv
        denom = s @ Binv_y

        if abs(denom) > 1e-30:
            u = s - Binv_y
            Binv += np.outer(u, sBinv) / denom

        x = x_new
        F_old = F_new
    else:
        iexit = 3

    return {
        'x': x,
        'fnorm': fnorm,
        'f': float(np.sum(F_old**2)),
        'niter': k if iexit != 0 else maxiter,
        'nfev': nfev,
        'iexit': iexit,
    }


def broyden_minimize(
    residual_func: Callable[[np.ndarray], np.ndarray],
    x0: np.ndarray,
    **kwargs,
) -> dict:
    """
    Broyden法で非線形最小二乗問題 min sum(r_i(x)^2) を解く。
    残差 r(x) = 0 を直接求解する（m = n の正方系が必要）。

    Parameters
    ----------
    residual_func : callable
        残差ベクトル r(x) を返す関数
    x0 : ndarray
        初期推定値
    **kwargs
        broyden_solve に渡す追加引数

    Returns
    -------
    dict
        broyden_solve の結果に 'f' (最終関数値 sum(r^2)) を追加
    """
    result = broyden_solve(residual_func, x0, **kwargs)
    r_final = residual_func(result['x'])
    result['f'] = float(np.sum(r_final**2))
    return result
