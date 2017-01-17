! ================================================================
! 
!> @note
!! This software was developed at the National Institute of
!! Standards and Technology by employees of the Federal Government
!! in the course of their official duties.  Pursuant to Title 17
!! Section 105 of the United States Code, this software is not
!! subject to copyright protection and is in the public domain.
!! This software is an experimental system.  NIST assumes no
!! responsibility whatsoever for its use by other parties, and makes
!! no guarantees, expressed or implied, about its quality,
!! reliability, or any other characteristic.  We would appreciate
!! acknowledgment if the software is used.
! 
! ================================================================
! 
!> @file
!! This file defines the math_utils module which groups some math
!! utility routines used by the rest of ZENO.
!!
!! @author Mark Mansfield (formerly of Stevens Institute) and students
!! @author Walid Keyrouz (walid.keyrouz@nist.gov)
!! @date   Mon Sep 23 14:16:36 2013 EDT
!! 
!! @defgroup math_utils Math Utilities
!! 
!! Groups together some math utility routines.
!!
!! @warning may be replaced by Fortran intrinsic functionality.
!! @warning routines to be renamed to be more descriptive of functionality.
!!
!! @todo check if these routines can be replaced by Fortran intrinsics.
! 
! Time-stamp: <2015-01-22 00:02:01 wtk>
! 
! ================================================================

module zeno_math_utils

  !! ================================================================

  use numeric_kinds

  implicit none
  private
  public summer, multiply, divide, power

contains

  ! ================================================================
  !
  ! Public routines
  !
  ! ================================================================

  !> @f$ c_1 = a_1 + b_1; c_2 = \sqrt{a_2^2 + b_2^2} @f$
  !! 
  !! @param a  (in)  double precision vector(2) 
  !! @param b  (in)  double precision vector(2)
  !! @param c  (out) double precision vector(2)

  subroutine summer(a, b, c)

    !! Declare arguments
    real(dp_k), dimension(2) :: a
    real(dp_k), dimension(2) :: b
    real(dp_k), dimension(2) :: c

    c(2) = a(2)**2 + b(2)**2
    c(2) = dsqrt(c(2))
    c(1) = a(1) + b(1)

  end subroutine summer

  !! ----------------------------------------------------------------

  !> @f$ c_1 = a_1 b_1; c_2 = \sqrt{b_1^2 a_2^2 + a_1^2 b_2^2} @f$
  !! 
  !! @param a  (in)  double precision vector, length = 2
  !! @param b  (in)  double precision vector, length = 2
  !! @param c  (out) double precision vector, length = 2

  subroutine multiply(a, b, c)

    !! Declare arguments

    real(dp_k), dimension(2) :: a, b, c
    
    c(2) = (b(1)**2)*(a(2)**2) + (a(1)**2)*(b(2)**2)
    c(2) = dsqrt(c(2))
    c(1) = a(1)*b(1)

  end subroutine multiply

  !! ----------------------------------------------------------------

  !> @brief to be filled in
  !! 
  !! @param x1 (in) double precision vector, length 2
  !! @param x2 (in) double precision vector, length 2
  !!
  !! @param y (out) double precision vector, length 2

  subroutine divide(x1, x2, y)

    !! Declare arguments
    real(dp_k), dimension(2) :: x1, x2, y

    !! Local variables
    real(dp_k) :: d1, d2
    
    d1 = 1.0d0/x2(1)
    d2 = -x1(1)/(x2(1)**2)

    y(2) = (d1**2)*(x1(2)**2) + (d2**2)*(x2(2)**2)
    y(2) = dsqrt(y(2))
    y(1) = x1(1)/x2(1)

  end subroutine divide

  !! ----------------------------------------------------------------

  !> @brief to be filled
  !! 
  !! @param x (in) double precision vector, length 2
  !! @param exponent (in) double precision scalar
  !!
  !! @param y (out) double precision vector, length 2

  subroutine power(x, y, exponent)

    !! Declare arguments
    real(dp_k), dimension(2) :: x, y
    real(dp_k)               :: exponent

    !! Local variables
    real(dp_k)               :: d
    
    y(1) = x(1)**exponent
    d = exponent * (x(1)**(exponent-1.0d0))

    y(2) = (d**2)*(x(2)**2)
    y(2) = dsqrt(y(2))

  end subroutine power

!! ================================================================

end module zeno_math_utils

!!! ================================================================

!!! Local Variables:
!!! mode: f90
!!! time-stamp-line-limit: 30
!!! fill-column: 72
!!! End:
