!========================================================================
! Levenberg-Marquardt法のテスト
! VA10A (minimize) と同一のテスト関数を残差形式で定義し比較する
!
! Test 1: Rosenbrock関数 (残差形式)
!   r1 = 1 - x,  r2 = 10*(y - x^2)
!   f = r1^2 + r2^2 = (1-x)^2 + 100*(y-x^2)^2
!   最小値: f(1, 1) = 0
!
! Test 2: 二次関数 (残差形式)
!   r1 = x - 3,  r2 = y + 2
!   f = r1^2 + r2^2 = (x-3)^2 + (y+2)^2
!   最小値: f(3, -2) = 0
!========================================================================
program test_lm
  use lm_solver_mod
  implicit none

  integer, parameter :: n = 2, m = 2
  real(dp) :: x(n), f
  integer :: iexit, npass, nfail
  real(dp) :: tol, err_x, err_y

  npass = 0
  nfail = 0

  write(6,*) '================================================'
  write(6,*) 'Levenberg-Marquardt solver test suite'
  write(6,*) '================================================'
  write(6,*)

  !--------------------------------------------------------------------
  ! Test 1: Rosenbrock
  !--------------------------------------------------------------------
  write(6,*) 'Test 1: Rosenbrock function (LM)'
  write(6,*) '  f(x,y) = (1-x)^2 + 100*(y-x^2)^2'
  write(6,*) '  Residuals: r1 = 1-x, r2 = 10*(y-x^2)'
  write(6,*) '  Expected minimum: f(1, 1) = 0'
  write(6,*)

  x(1) = -1.0_dp
  x(2) =  1.0_dp

  call lm_minimize(rosenbrock_resid, m, n, x, f, 1.0e-10_dp, 10000, 0, iexit)

  write(6,'(a,f12.6,a,f12.6)') '  Result: x = ', x(1), ', y = ', x(2)
  write(6,'(a,e15.7)') '  f(x,y) = ', f
  write(6,'(a,i3)') '  iexit  = ', iexit
  write(6,*)

  tol = 1.0e-4_dp
  err_x = abs(x(1) - 1.0_dp)
  err_y = abs(x(2) - 1.0_dp)
  if (err_x < tol .and. err_y < tol) then
    write(6,*) '  >>> PASS: Converged to correct minimum'
    npass = npass + 1
  else
    write(6,*) '  >>> FAIL: Did not converge to (1, 1)'
    write(6,'(a,e12.4,a,e12.4)') '  Error: |x-1| = ', err_x, ', |y-1| = ', err_y
    nfail = nfail + 1
  end if
  write(6,*)

  !--------------------------------------------------------------------
  ! Test 2: Quadratic
  !--------------------------------------------------------------------
  write(6,*) 'Test 2: Quadratic function (LM)'
  write(6,*) '  f(x,y) = (x-3)^2 + (y+2)^2'
  write(6,*) '  Residuals: r1 = x-3, r2 = y+2'
  write(6,*) '  Expected minimum: f(3, -2) = 0'
  write(6,*)

  x(1) = 0.0_dp
  x(2) = 0.0_dp

  call lm_minimize(quadratic_resid, m, n, x, f, 1.0e-10_dp, 10000, 0, iexit)

  write(6,'(a,f12.6,a,f12.6)') '  Result: x = ', x(1), ', y = ', x(2)
  write(6,'(a,e15.7)') '  f(x,y) = ', f
  write(6,'(a,i3)') '  iexit  = ', iexit
  write(6,*)

  tol = 1.0e-6_dp
  err_x = abs(x(1) - 3.0_dp)
  err_y = abs(x(2) + 2.0_dp)
  if (err_x < tol .and. err_y < tol) then
    write(6,*) '  >>> PASS: Converged to correct minimum'
    npass = npass + 1
  else
    write(6,*) '  >>> FAIL: Did not converge to (3, -2)'
    write(6,'(a,e12.4,a,e12.4)') '  Error: |x-3| = ', err_x, ', |y+2| = ', err_y
    nfail = nfail + 1
  end if
  write(6,*)

  !--------------------------------------------------------------------
  ! Summary
  !--------------------------------------------------------------------
  write(6,*) '================================================'
  write(6,'(a,i2,a,i2,a)') '  Results: ', npass, ' PASS, ', nfail, ' FAIL'
  write(6,*) '================================================'

  if (nfail > 0) stop 1

contains

  subroutine rosenbrock_resid(m, n, x, r)
    integer, intent(in) :: m, n
    real(dp), intent(in) :: x(n)
    real(dp), intent(out) :: r(m)
    r(1) = 1.0_dp - x(1)
    r(2) = 10.0_dp * (x(2) - x(1)**2)
  end subroutine

  subroutine quadratic_resid(m, n, x, r)
    integer, intent(in) :: m, n
    real(dp), intent(in) :: x(n)
    real(dp), intent(out) :: r(m)
    r(1) = x(1) - 3.0_dp
    r(2) = x(2) + 2.0_dp
  end subroutine

end program test_lm
