!========================================================================
! Broyden法（Broyden's first method / "good" Broyden）
! による非線形方程式系 F(x) = 0 のソルバー
!
! 最小化問題 min f(x) = sum r_i(x)^2 を、
! 非線形方程式 grad f(x) = 0 として解く。
! あるいは、残差 r(x) = 0 を直接解く（残差が正方の場合）。
!
! 参考文献:
!   C.G. Broyden, Math. Comp. 19, 577 (1965)
!   D.D. Johnson, Phys. Rev. B 38, 12807 (1988)
!========================================================================
module broyden_solver_mod
  implicit none
  integer, parameter :: dp = selected_real_kind(15, 307)

  abstract interface
    subroutine vector_func(n, x, f)
      import :: dp
      integer, intent(in) :: n
      real(dp), intent(in) :: x(n)
      real(dp), intent(out) :: f(n)
    end subroutine
  end interface

contains

  !----------------------------------------------------------------------
  ! Broyden法による非線形方程式 F(x) = 0 の求解
  !
  ! 引数:
  !   func    : F(x) を計算する関数
  !   n       : 変数/方程式の数
  !   x(n)    : 入出力。初期推定 → 解
  !   fnorm   : 出力。最終残差ノルム ||F(x)||
  !   eps     : 収束判定許容誤差
  !   maxiter : 最大反復回数
  !   iprint  : 出力制御
  !   iexit   : 出力。1=収束, 2=ステップ微小, 3=最大回数到達
  !----------------------------------------------------------------------
  subroutine broyden_solve(func, n, x, fnorm, eps, maxiter, iprint, iexit)
    procedure(vector_func) :: func
    integer, intent(in) :: n, maxiter, iprint
    real(dp), intent(inout) :: x(n)
    real(dp), intent(out) :: fnorm
    real(dp), intent(in) :: eps
    integer, intent(out) :: iexit

    real(dp) :: F_old(n), F_new(n), s(n), y(n), dx(n)
    real(dp) :: Binv(n, n)  ! Jacobian逆行列の近似
    real(dp) :: x_new(n)
    real(dp) :: hh, denom, step_norm
    integer :: iter, i, j

    ! 初期設定
    hh = 1.0e-7_dp
    iexit = 0

    ! 初期 Jacobian の数値計算と逆行列
    call init_jacobian_inv(func, n, x, hh, Binv)

    ! 初期残差
    call func(n, x, F_old)
    fnorm = sqrt(sum(F_old**2))

    if (iprint > 0) then
      write(6,'(a)') ' Entry into Broyden solver'
      write(6,'(a,i5,a,e15.7)') '  iter=', 0, '  ||F||=', fnorm
    end if

    do iter = 1, maxiter
      ! Newton-like ステップ: dx = -B^{-1} F
      do i = 1, n
        dx(i) = 0.0_dp
        do j = 1, n
          dx(i) = dx(i) - Binv(i, j) * F_old(j)
        end do
      end do

      ! ステップサイズの制限（安定化）
      step_norm = sqrt(sum(dx**2))
      if (step_norm > 10.0_dp) then
        dx = dx * (10.0_dp / step_norm)
        step_norm = 10.0_dp
      end if

      ! 更新
      x_new = x + dx
      call func(n, x_new, F_new)

      ! 残差ノルムの計算
      fnorm = sqrt(sum(F_new**2))

      if (iprint > 0) then
        if (mod(iter, iprint) == 0) then
          write(6,'(a,i5,a,e15.7)') '  iter=', iter, '  ||F||=', fnorm
        end if
      end if

      ! 収束判定
      if (fnorm < eps) then
        x = x_new
        iexit = 1
        exit
      end if

      if (step_norm < eps * (sqrt(sum(x**2)) + eps)) then
        x = x_new
        iexit = 2
        exit
      end if

      ! Broyden更新: B^{-1}_{k+1} = B^{-1}_k + (s - B^{-1}_k y) s^T B^{-1}_k / (s^T B^{-1}_k y)
      ! ここで s = dx, y = F_new - F_old
      s = dx
      y = F_new - F_old

      call update_broyden_inv(n, Binv, s, y)

      x = x_new
      F_old = F_new
    end do

    if (iexit == 0) iexit = 3

    if (iprint > 0) then
      write(6,'(a,i5,a,i3)') '  itn=', iter, '  iexit=', iexit
      write(6,'(a,e15.7)') '  ||F||=', fnorm
    end if

  end subroutine broyden_solve

  !----------------------------------------------------------------------
  ! 初期Jacobianの数値計算と逆行列（LU分解）
  !----------------------------------------------------------------------
  subroutine init_jacobian_inv(func, n, x, hh, Binv)
    procedure(vector_func) :: func
    integer, intent(in) :: n
    real(dp), intent(in) :: x(n), hh
    real(dp), intent(out) :: Binv(n, n)

    real(dp) :: F0(n), F1(n), x_pert(n)
    real(dp) :: Jac(n, n)
    integer :: jj

    call func(n, x, F0)

    do jj = 1, n
      x_pert = x
      x_pert(jj) = x(jj) + hh
      call func(n, x_pert, F1)
      Jac(:, jj) = (F1 - F0) / hh
    end do

    ! Jacobianの逆行列を計算（Gauss-Jordan）
    call invert_matrix(n, Jac, Binv)

  end subroutine init_jacobian_inv

  !----------------------------------------------------------------------
  ! Broyden逆行列更新（Sherman-Morrison公式）
  ! B^{-1}_{k+1} = B^{-1}_k + (s - B^{-1}_k y) * (s^T B^{-1}_k) / (s^T B^{-1}_k y)
  !----------------------------------------------------------------------
  subroutine update_broyden_inv(n, Binv, s, y)
    integer, intent(in) :: n
    real(dp), intent(inout) :: Binv(n, n)
    real(dp), intent(in) :: s(n), y(n)

    real(dp) :: Binv_y(n), sBinv(n), u(n)
    real(dp) :: denom
    integer :: i, j

    ! Binv * y
    do i = 1, n
      Binv_y(i) = dot_product(Binv(i, :), y)
    end do

    ! s^T * Binv
    do j = 1, n
      sBinv(j) = dot_product(s, Binv(:, j))
    end do

    ! s^T * Binv * y
    denom = dot_product(s, Binv_y)

    if (abs(denom) < 1.0e-30_dp) return  ! 特異な場合はスキップ

    ! u = s - Binv * y
    u = s - Binv_y

    ! Binv += u * sBinv^T / denom
    do i = 1, n
      do j = 1, n
        Binv(i, j) = Binv(i, j) + u(i) * sBinv(j) / denom
      end do
    end do

  end subroutine update_broyden_inv

  !----------------------------------------------------------------------
  ! 行列の逆行列（Gauss-Jordan法）
  !----------------------------------------------------------------------
  subroutine invert_matrix(n, A, Ainv)
    integer, intent(in) :: n
    real(dp), intent(in) :: A(n, n)
    real(dp), intent(out) :: Ainv(n, n)

    real(dp) :: Aug(n, 2*n), pivot, factor
    integer :: i, j, k, max_row

    ! 拡大行列 [A | I] の構成
    Aug(:, 1:n) = A
    Aug(:, n+1:2*n) = 0.0_dp
    do i = 1, n
      Aug(i, n+i) = 1.0_dp
    end do

    ! 前進消去（部分ピボット選択付き）
    do k = 1, n
      ! ピボット選択
      max_row = k
      do i = k + 1, n
        if (abs(Aug(i, k)) > abs(Aug(max_row, k))) max_row = i
      end do
      if (max_row /= k) then
        do j = 1, 2*n
          factor = Aug(k, j)
          Aug(k, j) = Aug(max_row, j)
          Aug(max_row, j) = factor
        end do
      end if

      pivot = Aug(k, k)
      if (abs(pivot) < 1.0e-30_dp) then
        ! 特異行列：単位行列で代用
        Ainv = 0.0_dp
        do i = 1, n
          Ainv(i, i) = 1.0_dp
        end do
        return
      end if

      Aug(k, :) = Aug(k, :) / pivot

      do i = 1, n
        if (i /= k) then
          factor = Aug(i, k)
          Aug(i, :) = Aug(i, :) - factor * Aug(k, :)
        end if
      end do
    end do

    Ainv = Aug(:, n+1:2*n)
  end subroutine invert_matrix

  !----------------------------------------------------------------------
  ! 最小化問題のラッパー：
  ! min f(x) = sum r_i(x)^2 を F(x) = 2 J^T r = 0 として解く代わりに、
  ! 残差 r(x) = 0 を直接解く（m = n の場合）
  !----------------------------------------------------------------------
  subroutine broyden_minimize(resid, m, n, x, f, eps, maxiter, iprint, iexit)
    procedure(vector_func) :: resid
    integer, intent(in) :: m, n, maxiter, iprint
    real(dp), intent(inout) :: x(n)
    real(dp), intent(out) :: f
    real(dp), intent(in) :: eps
    integer, intent(out) :: iexit

    real(dp) :: fnorm, r(n)

    ! m = n が必要（正方系）
    if (m /= n) then
      write(6,*) 'ERROR: broyden_minimize requires m == n'
      iexit = 3
      f = huge(1.0_dp)
      return
    end if

    call broyden_solve(resid, n, x, fnorm, eps, maxiter, iprint, iexit)

    ! 最終関数値の計算
    call resid(n, x, r)
    f = sum(r**2)

  end subroutine broyden_minimize

end module broyden_solver_mod
