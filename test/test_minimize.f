C========+=========+=========+=========+=========+=========+=========+=$
C PROGRAM: test_minimize
C PURPOSE: Test the minimize subroutine with two nonlinear functions
C          with analytically known solutions.
C
C Test 1: Rosenbrock function  f(x,y) = (1-x)^2 + 100*(y-x^2)^2
C         Minimum at (1, 1) with f = 0
C
C Test 2: Quadratic function   f(x,y) = (x-3)^2 + (y+2)^2
C         Minimum at (3, -2) with f = 0
C========+=========+=========+=========+=========+=========+=========+=$
      program test_minimize
      implicit double precision (a-h,o-z)
      parameter (nmax=10, nhmax=nmax*(nmax+1)/2)
      parameter (nwmax=nmax*nmax)
      dimension x(nmax), g(nmax), h(nhmax), w(nwmax), xm(nmax)
      external rosenbrock, quadratic

      integer npass, nfail
      npass = 0
      nfail = 0

      write(6,*) '================================================'
      write(6,*) 'minimize subroutine test suite'
      write(6,*) '================================================'
      write(6,*)

C----------------------------------------------------------------------
C Test 1: Rosenbrock function
C----------------------------------------------------------------------
      write(6,*) 'Test 1: Rosenbrock function'
      write(6,*) '  f(x,y) = (1-x)^2 + 100*(y-x^2)^2'
      write(6,*) '  Expected minimum: f(1, 1) = 0'
      write(6,*)

      n = 2
      x(1) = -1.0d0
      x(2) =  1.0d0
      f = 0.0d0

      do i = 1, n
         xm(i) = dabs(x(i)) + 1.0d-15
      end do

      mode = 1
      dfn = -0.5d0
      hh = 1.0d-5
      eps = 1.0d-5
      maxfn = 10000
      iprint = 0
      iexit = 0

      call minimize(rosenbrock, n, x, f, g, h, w, dfn, xm,
     $              hh, eps, mode, maxfn, iprint, iexit)

      write(6,100) x(1), x(2)
 100  format('  Result: x = ', f12.6, ', y = ', f12.6)
      write(6,101) f
 101  format('  f(x,y) = ', e15.7)
      write(6,102) iexit
 102  format('  iexit  = ', i3)
      write(6,*)

      tol = 1.0d-2
      err_x = dabs(x(1) - 1.0d0)
      err_y = dabs(x(2) - 1.0d0)
      if (err_x .lt. tol .and. err_y .lt. tol) then
         write(6,*) '  >>> PASS: Converged to correct minimum'
         npass = npass + 1
      else
         write(6,*) '  >>> FAIL: Did not converge to (1, 1)'
         write(6,103) err_x, err_y
 103     format('  Error: |x-1| = ', e12.4, ', |y-1| = ', e12.4)
         nfail = nfail + 1
      end if
      write(6,*)

C----------------------------------------------------------------------
C Test 2: Quadratic function
C----------------------------------------------------------------------
      write(6,*) 'Test 2: Quadratic function'
      write(6,*) '  f(x,y) = (x-3)^2 + (y+2)^2'
      write(6,*) '  Expected minimum: f(3, -2) = 0'
      write(6,*)

      n = 2
      x(1) = 0.0d0
      x(2) = 0.0d0
      f = 0.0d0

      do i = 1, n
         xm(i) = dabs(x(i)) + 1.0d0
      end do

      mode = 1
      dfn = -0.5d0
      hh = 1.0d-5
      eps = 1.0d-5
      maxfn = 10000
      iprint = 0
      iexit = 0

      call minimize(quadratic, n, x, f, g, h, w, dfn, xm,
     $              hh, eps, mode, maxfn, iprint, iexit)

      write(6,200) x(1), x(2)
 200  format('  Result: x = ', f12.6, ', y = ', f12.6)
      write(6,201) f
 201  format('  f(x,y) = ', e15.7)
      write(6,202) iexit
 202  format('  iexit  = ', i3)
      write(6,*)

      tol = 1.0d-3
      err_x = dabs(x(1) - 3.0d0)
      err_y = dabs(x(2) + 2.0d0)
      if (err_x .lt. tol .and. err_y .lt. tol) then
         write(6,*) '  >>> PASS: Converged to correct minimum'
         npass = npass + 1
      else
         write(6,*) '  >>> FAIL: Did not converge to (3, -2)'
         write(6,203) err_x, err_y
 203     format('  Error: |x-3| = ', e12.4, ', |y+2| = ', e12.4)
         nfail = nfail + 1
      end if
      write(6,*)

C----------------------------------------------------------------------
C Summary
C----------------------------------------------------------------------
      write(6,*) '================================================'
      write(6,300) npass, nfail
 300  format(' Results: ', i2, ' PASS, ', i2, ' FAIL')
      write(6,*) '================================================'

      if (nfail .gt. 0) then
         stop 1
      end if

      end

C========+=========+=========+=========+=========+=========+=========+=$
C Rosenbrock function: f(x,y) = (1-x)^2 + 100*(y-x^2)^2
C Minimum at (1, 1), f = 0
C========+=========+=========+=========+=========+=========+=========+=$
      subroutine rosenbrock(n, x, f)
      implicit double precision (a-h,o-z)
      dimension x(*)
      f = (1.0d0 - x(1))**2 + 100.0d0*(x(2) - x(1)**2)**2
      return
      end

C========+=========+=========+=========+=========+=========+=========+=$
C Quadratic function: f(x,y) = (x-3)^2 + (y+2)^2
C Minimum at (3, -2), f = 0
C========+=========+=========+=========+=========+=========+=========+=$
      subroutine quadratic(n, x, f)
      implicit double precision (a-h,o-z)
      dimension x(*)
      f = (x(1) - 3.0d0)**2 + (x(2) + 2.0d0)**2
      return
      end
