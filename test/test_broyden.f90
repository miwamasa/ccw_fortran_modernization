!========================================================================
! Broyden法のテスト
! 残差 r(x) = 0 を直接解く（m = n の正方系）
!
! Test 1: Rosenbrock関数
!   r1 = 1 - x = 0  →  x = 1
!   r2 = 10*(y - x^2) = 0  →  y = x^2 = 1
!   f = r1^2 + r2^2 = (1-x)^2 + 100*(y-x^2)^2
!
! Test 2: 二次関数
!   r1 = x - 3 = 0  →  x = 3
!   r2 = y + 2 = 0  →  y = -2
!   f = r1^2 + r2^2 = (x-3)^2 + (y+2)^2
!========================================================================
program test_broyden
  use broyden_solver_mod
  implicit none

  integer, parameter :: n = 2
  real(dp) :: x(n), f, fnorm
  integer :: iexit, npass, nfail
  real(dp) :: tol, err_x, err_y

  npass = 0
  nfail = 0

  write(6,*) '================================================'
  write(6,*) 'Broyden solver test suite'
  write(6,*) '================================================'
  write(6,*)

  !--------------------------------------------------------------------
  ! Test 1: Rosenbrock
  !--------------------------------------------------------------------
  write(6,*) 'Test 1: Rosenbrock function (Broyden)'
  write(6,*) '  r1 = 1-x, r2 = 10*(y-x^2)'
  write(6,*) '  Solving r(x) = 0 directly'
  write(6,*) '  Expected solution: (1, 1)'
  write(6,*)

  x(1) = -1.0_dp
  x(2) =  1.0_dp

  call broyden_minimize(rosenbrock_resid, n, n, x, f, 1.0e-10_dp, 10000, 0, iexit)

  write(6,'(a,f12.6,a,f12.6)') '  Result: x = ', x(1), ', y = ', x(2)
  write(6,'(a,e15.7)') '  f(x,y) = ', f
  write(6,'(a,i3)') '  iexit  = ', iexit
  write(6,*)

  tol = 1.0e-4_dp
  err_x = abs(x(1) - 1.0_dp)
  err_y = abs(x(2) - 1.0_dp)
  if (err_x < tol .and. err_y < tol) then
    write(6,*) '  >>> PASS: Converged to correct solution'
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
  write(6,*) 'Test 2: Quadratic function (Broyden)'
  write(6,*) '  r1 = x-3, r2 = y+2'
  write(6,*) '  Solving r(x) = 0 directly'
  write(6,*) '  Expected solution: (3, -2)'
  write(6,*)

  x(1) = 0.0_dp
  x(2) = 0.0_dp

  call broyden_minimize(quadratic_resid, n, n, x, f, 1.0e-10_dp, 10000, 0, iexit)

  write(6,'(a,f12.6,a,f12.6)') '  Result: x = ', x(1), ', y = ', x(2)
  write(6,'(a,e15.7)') '  f(x,y) = ', f
  write(6,'(a,i3)') '  iexit  = ', iexit
  write(6,*)

  tol = 1.0e-6_dp
  err_x = abs(x(1) - 3.0_dp)
  err_y = abs(x(2) + 2.0_dp)
  if (err_x < tol .and. err_y < tol) then
    write(6,*) '  >>> PASS: Converged to correct solution'
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

  subroutine rosenbrock_resid(n, x, r)
    integer, intent(in) :: n
    real(dp), intent(in) :: x(n)
    real(dp), intent(out) :: r(n)
    r(1) = 1.0_dp - x(1)
    r(2) = 10.0_dp * (x(2) - x(1)**2)
  end subroutine

  subroutine quadratic_resid(n, x, r)
    integer, intent(in) :: n
    real(dp), intent(in) :: x(n)
    real(dp), intent(out) :: r(n)
    r(1) = x(1) - 3.0_dp
    r(2) = x(2) + 2.0_dp
  end subroutine

end program test_broyden
