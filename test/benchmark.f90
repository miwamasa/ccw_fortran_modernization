!========================================================================
! 性能比較ベンチマーク: VA10A vs Levenberg-Marquardt vs Broyden
!
! テスト関数:
!   1. Rosenbrock 2D  (簡単)
!   2. Extended Rosenbrock 6D (中程度)
!   3. Trigonometric 10D (困難)
!
! 各関数について反復数、関数評価回数、精度、実行時間を計測
!========================================================================
program benchmark
  use lm_solver_mod, only: dp, lm_minimize, residual_func
  use broyden_solver_mod, only: broyden_minimize, vector_func
  implicit none

  integer, parameter :: nrep = 100  ! 時間計測の繰り返し回数
  real(dp) :: t_start, t_end

  ! 外部関数の宣言
  interface
    subroutine rosenbrock_lm(m, n, x, r)
      import :: dp
      integer, intent(in) :: m, n
      real(dp), intent(in) :: x(n)
      real(dp), intent(out) :: r(m)
    end subroutine
    subroutine rosenbrock6_lm(m, n, x, r)
      import :: dp
      integer, intent(in) :: m, n
      real(dp), intent(in) :: x(n)
      real(dp), intent(out) :: r(m)
    end subroutine
    subroutine trig10_lm(m, n, x, r)
      import :: dp
      integer, intent(in) :: m, n
      real(dp), intent(in) :: x(n)
      real(dp), intent(out) :: r(m)
    end subroutine
    subroutine rosenbrock_br(n, x, r)
      import :: dp
      integer, intent(in) :: n
      real(dp), intent(in) :: x(n)
      real(dp), intent(out) :: r(n)
    end subroutine
    subroutine trig10_br(n, x, r)
      import :: dp
      integer, intent(in) :: n
      real(dp), intent(in) :: x(n)
      real(dp), intent(out) :: r(n)
    end subroutine
  end interface

  write(6,'(a)') '================================================================'
  write(6,'(a)') ' Performance Benchmark: VA10A vs Levenberg-Marquardt vs Broyden'
  write(6,'(a)') '================================================================'
  write(6,*)

  call bench_rosenbrock_2d()
  call bench_rosenbrock_6d()
  call bench_trig_10d()

contains

  !--------------------------------------------------------------------
  ! Benchmark 1: Rosenbrock 2D
  !--------------------------------------------------------------------
  subroutine bench_rosenbrock_2d()
    implicit none
    integer, parameter :: n = 2
    real(dp) :: x(n), f, g(n), h(n*(n+1)/2), w(n*n*4), xm(n)
    real(dp) :: dfn, hh_va, deps
    integer :: mode, maxfn, iprint_va, iexit, i, rep
    external :: rosenbrock_va

    write(6,'(a)') '--- Benchmark 1: Rosenbrock 2D ---'
    write(6,'(a)') '    f(x,y) = (1-x)^2 + 100*(y-x^2)^2'
    write(6,'(a)') '    Initial: (-1, 1), Solution: (1, 1)'
    write(6,*)

    ! --- VA10A ---
    call cpu_time(t_start)
    do rep = 1, nrep
      x = (/ -1.0_dp, 1.0_dp /)
      f = 0.0_dp
      do i = 1, n
        xm(i) = abs(x(i)) + 1.0e-15_dp
      end do
      mode = 1; dfn = -0.5_dp; hh_va = 1.0e-5_dp; deps = 1.0e-5_dp
      maxfn = 10000; iprint_va = 0; iexit = 0
      call minimize(rosenbrock_va, n, x, f, g, h, w, dfn, xm, &
                    hh_va, deps, mode, maxfn, iprint_va, iexit)
    end do
    call cpu_time(t_end)
    write(6,'(a,f12.6,a,f12.6,a,e12.4,a,i2,a,f8.4,a)') &
      '  VA10A:  x=', x(1), ' y=', x(2), ' f=', f, ' iexit=', iexit, &
      ' time=', (t_end-t_start)/nrep*1000, 'ms'

    ! --- Levenberg-Marquardt ---
    call cpu_time(t_start)
    do rep = 1, nrep
      x = (/ -1.0_dp, 1.0_dp /)
      call lm_minimize(rosenbrock_lm, 2, 2, x, f, 1.0e-10_dp, 10000, 0, iexit)
    end do
    call cpu_time(t_end)
    write(6,'(a,f12.6,a,f12.6,a,e12.4,a,i2,a,f8.4,a)') &
      '  LM:     x=', x(1), ' y=', x(2), ' f=', f, ' iexit=', iexit, &
      ' time=', (t_end-t_start)/nrep*1000, 'ms'

    ! --- Broyden ---
    call cpu_time(t_start)
    do rep = 1, nrep
      x = (/ -1.0_dp, 1.0_dp /)
      call broyden_minimize(rosenbrock_br, 2, 2, x, f, 1.0e-10_dp, 10000, 0, iexit)
    end do
    call cpu_time(t_end)
    write(6,'(a,f12.6,a,f12.6,a,e12.4,a,i2,a,f8.4,a)') &
      '  Broyden:x=', x(1), ' y=', x(2), ' f=', f, ' iexit=', iexit, &
      ' time=', (t_end-t_start)/nrep*1000, 'ms'

    write(6,*)
  end subroutine

  !--------------------------------------------------------------------
  ! Benchmark 2: Extended Rosenbrock 6D
  !--------------------------------------------------------------------
  subroutine bench_rosenbrock_6d()
    implicit none
    integer, parameter :: n = 6
    real(dp) :: x(n), f, g(n), h(n*(n+1)/2), w(n*n*4), xm(n)
    real(dp) :: dfn, hh_va, deps
    integer :: mode, maxfn, iprint_va, iexit, i, rep
    external :: rosenbrock6_va

    write(6,'(a)') '--- Benchmark 2: Extended Rosenbrock 6D ---'
    write(6,'(a)') '    f(x) = sum_{i=1}^{5} [(1-x_i)^2 + 100*(x_{i+1}-x_i^2)^2]'
    write(6,'(a)') '    Initial: (-1.2,1,-1.2,1,-1.2,1), Solution: (1,...,1)'
    write(6,*)

    ! --- VA10A ---
    call cpu_time(t_start)
    do rep = 1, nrep
      x = (/ -1.2_dp, 1.0_dp, -1.2_dp, 1.0_dp, -1.2_dp, 1.0_dp /)
      f = 0.0_dp
      do i = 1, n
        xm(i) = abs(x(i)) + 1.0e-15_dp
      end do
      mode = 1; dfn = -0.5_dp; hh_va = 1.0e-5_dp; deps = 1.0e-5_dp
      maxfn = 50000; iprint_va = 0; iexit = 0
      call minimize(rosenbrock6_va, n, x, f, g, h, w, dfn, xm, &
                    hh_va, deps, mode, maxfn, iprint_va, iexit)
    end do
    call cpu_time(t_end)
    write(6,'(a,e12.4,a,e12.4,a,i2,a,f8.4,a)') &
      '  VA10A:  |x-1|_max=', maxval(abs(x-1.0_dp)), ' f=', f, &
      ' iexit=', iexit, ' time=', (t_end-t_start)/nrep*1000, 'ms'

    ! --- LM ---
    call cpu_time(t_start)
    do rep = 1, nrep
      x = (/ -1.2_dp, 1.0_dp, -1.2_dp, 1.0_dp, -1.2_dp, 1.0_dp /)
      call lm_minimize(rosenbrock6_lm, 5, 6, x, f, 1.0e-10_dp, 50000, 0, iexit)
    end do
    call cpu_time(t_end)
    write(6,'(a,e12.4,a,e12.4,a,i2,a,f8.4,a)') &
      '  LM:     |x-1|_max=', maxval(abs(x-1.0_dp)), ' f=', f, &
      ' iexit=', iexit, ' time=', (t_end-t_start)/nrep*1000, 'ms'

    write(6,*)
  end subroutine

  !--------------------------------------------------------------------
  ! Benchmark 3: Trigonometric 10D
  ! f(x) = sum_{i=1}^{10} [n - sum_j cos(x_j) + i*(1-cos(x_i)) - sin(x_i)]^2
  !--------------------------------------------------------------------
  subroutine bench_trig_10d()
    implicit none
    integer, parameter :: n = 10
    real(dp) :: x(n), f, g(n), h(n*(n+1)/2), w(n*n*4), xm(n)
    real(dp) :: dfn, hh_va, deps
    integer :: mode, maxfn, iprint_va, iexit, i, rep

    external :: trig10_va

    write(6,'(a)') '--- Benchmark 3: Trigonometric 10D ---'
    write(6,'(a)') '    r_i = n - sum_j cos(x_j) + i*(1-cos(x_i)) - sin(x_i)'
    write(6,'(a)') '    Initial: (1/n,...,1/n)'
    write(6,*)

    ! --- VA10A ---
    call cpu_time(t_start)
    do rep = 1, nrep
      do i = 1, n
        x(i) = 1.0_dp / n
      end do
      f = 0.0_dp
      do i = 1, n
        xm(i) = abs(x(i)) + 1.0e-15_dp
      end do
      mode = 1; dfn = -0.5_dp; hh_va = 1.0e-5_dp; deps = 1.0e-6_dp
      maxfn = 50000; iprint_va = 0; iexit = 0
      call minimize(trig10_va, n, x, f, g, h, w, dfn, xm, &
                    hh_va, deps, mode, maxfn, iprint_va, iexit)
    end do
    call cpu_time(t_end)
    write(6,'(a,e12.4,a,i2,a,f8.4,a)') &
      '  VA10A:  f=', f, ' iexit=', iexit, &
      ' time=', (t_end-t_start)/nrep*1000, 'ms'

    ! --- LM ---
    call cpu_time(t_start)
    do rep = 1, nrep
      do i = 1, n
        x(i) = 1.0_dp / n
      end do
      call lm_minimize(trig10_lm, 10, 10, x, f, 1.0e-12_dp, 50000, 0, iexit)
    end do
    call cpu_time(t_end)
    write(6,'(a,e12.4,a,i2,a,f8.4,a)') &
      '  LM:     f=', f, ' iexit=', iexit, &
      ' time=', (t_end-t_start)/nrep*1000, 'ms'

    ! --- Broyden ---
    call cpu_time(t_start)
    do rep = 1, nrep
      do i = 1, n
        x(i) = 1.0_dp / n
      end do
      call broyden_minimize(trig10_br, 10, 10, x, f, 1.0e-10_dp, 50000, 0, iexit)
    end do
    call cpu_time(t_end)
    write(6,'(a,e12.4,a,i2,a,f8.4,a)') &
      '  Broyden:f=', f, ' iexit=', iexit, &
      ' time=', (t_end-t_start)/nrep*1000, 'ms'

    write(6,*)
  end subroutine

end program benchmark

!======================================================================
! VA10A用の目的関数（subroutine(n, x, f)形式）
!======================================================================
subroutine rosenbrock_va(n, x, f)
  implicit double precision (a-h,o-z)
  dimension x(*)
  f = (1.0d0 - x(1))**2 + 100.0d0*(x(2) - x(1)**2)**2
end subroutine

subroutine rosenbrock6_va(n, x, f)
  implicit double precision (a-h,o-z)
  dimension x(*)
  f = 0.0d0
  do i = 1, n-1
    f = f + (1.0d0-x(i))**2 + 100.0d0*(x(i+1)-x(i)**2)**2
  end do
end subroutine

subroutine trig10_va(n, x, f)
  implicit double precision (a-h,o-z)
  dimension x(*)
  f = 0.0d0
  csum = 0.0d0
  do j = 1, n
    csum = csum + cos(x(j))
  end do
  do i = 1, n
    ri = dble(n) - csum + dble(i)*(1.0d0 - cos(x(i))) - sin(x(i))
    f = f + ri*ri
  end do
end subroutine

!======================================================================
! LM用の残差関数（subroutine(m, n, x, r)形式）
!======================================================================
subroutine rosenbrock_lm(m, n, x, r)
  use lm_solver_mod, only: dp
  integer, intent(in) :: m, n
  real(dp), intent(in) :: x(n)
  real(dp), intent(out) :: r(m)
  r(1) = 1.0_dp - x(1)
  r(2) = 10.0_dp * (x(2) - x(1)**2)
end subroutine

subroutine rosenbrock6_lm(m, n, x, r)
  use lm_solver_mod, only: dp
  integer, intent(in) :: m, n
  real(dp), intent(in) :: x(n)
  real(dp), intent(out) :: r(m)
  integer :: i
  do i = 1, m
    r(i) = 0.0_dp
  end do
  ! Extended Rosenbrock: r_{2i-1} = 1-x_i, r_{2i} = 10*(x_{i+1}-x_i^2)
  ! But m=5, n=6, so: r_i = (1-x_i)^2 + 100*(x_{i+1}-x_i^2)^2 won't work for LM directly
  ! Use: r_i = 10*(x_{i+1} - x_i^2)  and handle (1-x_i) separately
  ! Actually for m=5, n=6: each r_i combines both terms
  ! Let's use sqrt of each pair:
  do i = 1, 5
    r(i) = sqrt((1.0_dp-x(i))**2 + 100.0_dp*(x(i+1)-x(i)**2)**2)
  end do
end subroutine

subroutine trig10_lm(m, n, x, r)
  use lm_solver_mod, only: dp
  integer, intent(in) :: m, n
  real(dp), intent(in) :: x(n)
  real(dp), intent(out) :: r(m)
  real(dp) :: csum
  integer :: i, j
  csum = 0.0_dp
  do j = 1, n
    csum = csum + cos(x(j))
  end do
  do i = 1, m
    r(i) = dble(n) - csum + dble(i)*(1.0_dp - cos(x(i))) - sin(x(i))
  end do
end subroutine

!======================================================================
! Broyden用の残差関数（subroutine(n, x, r)形式）
!======================================================================
subroutine rosenbrock_br(n, x, r)
  use broyden_solver_mod, only: dp
  integer, intent(in) :: n
  real(dp), intent(in) :: x(n)
  real(dp), intent(out) :: r(n)
  r(1) = 1.0_dp - x(1)
  r(2) = 10.0_dp * (x(2) - x(1)**2)
end subroutine

subroutine trig10_br(n, x, r)
  use broyden_solver_mod, only: dp
  integer, intent(in) :: n
  real(dp), intent(in) :: x(n)
  real(dp), intent(out) :: r(n)
  real(dp) :: csum
  integer :: i, j
  csum = 0.0_dp
  do j = 1, n
    csum = csum + cos(x(j))
  end do
  do i = 1, n
    r(i) = dble(n) - csum + dble(i)*(1.0_dp - cos(x(i))) - sin(x(i))
  end do
end subroutine
