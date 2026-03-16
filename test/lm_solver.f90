!========================================================================
! Levenberg-Marquardt法による非線形最小二乗ソルバー
!
! 目的関数 f(x) = sum_i r_i(x)^2 を最小化する。
! Jacobian は数値微分で計算。
!
! 参考文献:
!   K. Levenberg, Quart. Appl. Math. 2, 164 (1944)
!   D.W. Marquardt, SIAM J. Appl. Math. 11, 431 (1963)
!========================================================================
module lm_solver_mod
  implicit none
  integer, parameter :: dp = selected_real_kind(15, 307)

  abstract interface
    subroutine residual_func(m, n, x, r)
      import :: dp
      integer, intent(in) :: m, n
      real(dp), intent(in) :: x(n)
      real(dp), intent(out) :: r(m)
    end subroutine
  end interface

contains

  subroutine lm_minimize(resid, m, n, x, f, eps, maxiter, iprint, iexit)
    procedure(residual_func) :: resid
    integer, intent(in) :: m, n, maxiter, iprint
    real(dp), intent(inout) :: x(n)
    real(dp), intent(out) :: f
    real(dp), intent(in) :: eps
    integer, intent(out) :: iexit

    real(dp) :: r(m), r_trial(m), Jac(m, n), x_trial(n)
    real(dp) :: A(n, n), g(n), delta(n), b(n)
    real(dp) :: lambda, f_trial, hh, rho, actual, predicted
    real(dp) :: pivot, factor, tmp
    integer :: iter, i, j, k, ifn, max_row, ii

    hh = 1.0e-7_dp
    lambda = 1.0e-3_dp
    iexit = 0
    ifn = 0

    ! 初期残差
    call resid(m, n, x, r)
    ifn = ifn + 1
    f = 0.0_dp
    do i = 1, m
      f = f + r(i) * r(i)
    end do

    if (iprint > 0) then
      write(6,'(a)') ' Entry into Levenberg-Marquardt solver'
      write(6,'(a,i5,a,e15.7)') '  iter=', 0, '  f=', f
    end if

    do iter = 1, maxiter
      ! Jacobian の数値計算（前進差分）
      do j = 1, n
        x_trial = x
        x_trial(j) = x_trial(j) + hh
        call resid(m, n, x_trial, r_trial)
        ifn = ifn + 1
        do i = 1, m
          Jac(i, j) = (r_trial(i) - r(i)) / hh
        end do
      end do

      ! A = J^T J,  g = J^T r
      do i = 1, n
        do j = 1, n
          A(i, j) = 0.0_dp
          do k = 1, m
            A(i, j) = A(i, j) + Jac(k, i) * Jac(k, j)
          end do
        end do
        g(i) = 0.0_dp
        do k = 1, m
          g(i) = g(i) + Jac(k, i) * r(k)
        end do
      end do

      ! 初回のlambdaスケーリング
      if (iter == 1) then
        tmp = 0.0_dp
        do i = 1, n
          if (A(i, i) > tmp) tmp = A(i, i)
        end do
        if (tmp < 1.0e-20_dp) tmp = 1.0_dp
        lambda = 1.0e-3_dp * tmp
      end if

      ! 勾配チェック（収束判定）
      tmp = 0.0_dp
      do i = 1, n
        if (abs(g(i)) > tmp) tmp = abs(g(i))
      end do
      if (tmp < eps) then
        iexit = 1
        exit
      end if

      ! LMステップ試行（lambda調整ループ）
      do k = 1, 30
        ! (J^T J + lambda * I) delta = -J^T r をGauss消去で解く
        b = -g
        do i = 1, n
          do j = 1, n
            x_trial(1) = A(i, j)  ! temp use
          end do
        end do
        ! Aのコピーを作って解く
        call solve_damped_system(n, A, g, lambda, delta)

        ! 試行点
        x_trial = x + delta
        call resid(m, n, x_trial, r_trial)
        ifn = ifn + 1
        f_trial = 0.0_dp
        do i = 1, m
          f_trial = f_trial + r_trial(i) * r_trial(i)
        end do

        ! 実際の改善量
        actual = f - f_trial

        ! 予測改善量: -delta^T (g + 0.5 * A * delta) ≈ -delta^T g - 0.5 lambda ||delta||^2
        predicted = 0.0_dp
        do i = 1, n
          predicted = predicted - delta(i) * g(i) + 0.5_dp * lambda * delta(i) * delta(i)
        end do

        if (predicted < 1.0e-30_dp) predicted = 1.0e-30_dp
        rho = actual / predicted

        if (rho > 1.0e-4_dp) then
          ! ステップ採用
          x = x_trial
          r = r_trial
          f = f_trial
          lambda = lambda * max(0.333_dp, 1.0_dp - (2.0_dp * rho - 1.0_dp)**3)
          if (lambda < 1.0e-20_dp) lambda = 1.0e-20_dp
          exit
        else
          ! ステップ棄却、lambdaを増大
          lambda = lambda * 4.0_dp
          if (lambda > 1.0e20_dp) lambda = 1.0e20_dp
        end if
      end do

      if (iprint > 0) then
        if (mod(iter, iprint) == 0) then
          write(6,'(a,i5,a,i5,a,e15.7,a,e10.3)') '  iter=', iter, &
            '  ifn=', ifn, '  f=', f, '  lam=', lambda
        end if
      end if

      ! 関数値による収束判定
      if (f < eps * eps) then
        iexit = 1
        exit
      end if

      ! ステップサイズによる収束判定
      tmp = 0.0_dp
      do i = 1, n
        if (abs(delta(i)) > tmp) tmp = abs(delta(i))
      end do
      if (tmp < eps * (maxval(abs(x)) + eps)) then
        iexit = 2
        exit
      end if
    end do

    if (iexit == 0) iexit = 3

    if (iprint > 0) then
      write(6,'(a,i5,a,i5,a,i3)') '  itn=', min(iter, maxiter), &
        '  ifn=', ifn, '  iexit=', iexit
      write(6,'(a,e15.7)') '  f=', f
    end if

  end subroutine lm_minimize

  !----------------------------------------------------------------------
  ! (A + lambda * I) delta = -g を部分ピボットGauss消去で解く
  !----------------------------------------------------------------------
  subroutine solve_damped_system(n, A_in, g, lambda, delta)
    integer, intent(in) :: n
    real(dp), intent(in) :: A_in(n, n), g(n), lambda
    real(dp), intent(out) :: delta(n)

    real(dp) :: A(n, n), b(n), pivot, factor, tmp
    integer :: i, k, max_row, ii

    ! Aのコピーにダンピング追加
    A = A_in
    do i = 1, n
      A(i, i) = A(i, i) + lambda
    end do
    b = -g

    ! 前進消去（部分ピボット選択）
    do k = 1, n
      max_row = k
      do ii = k + 1, n
        if (abs(A(ii, k)) > abs(A(max_row, k))) max_row = ii
      end do
      if (max_row /= k) then
        do ii = 1, n
          tmp = A(k, ii); A(k, ii) = A(max_row, ii); A(max_row, ii) = tmp
        end do
        tmp = b(k); b(k) = b(max_row); b(max_row) = tmp
      end if

      pivot = A(k, k)
      if (abs(pivot) < 1.0e-30_dp) then
        delta = 0.0_dp
        return
      end if

      do i = k + 1, n
        factor = A(i, k) / pivot
        do ii = k + 1, n
          A(i, ii) = A(i, ii) - factor * A(k, ii)
        end do
        b(i) = b(i) - factor * b(k)
      end do
    end do

    ! 後退代入
    do i = n, 1, -1
      delta(i) = b(i)
      do k = i + 1, n
        delta(i) = delta(i) - A(i, k) * delta(k)
      end do
      delta(i) = delta(i) / A(i, i)
    end do
  end subroutine solve_damped_system

end module lm_solver_mod
