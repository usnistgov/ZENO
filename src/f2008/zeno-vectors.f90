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
!! This file provides the `zeno_vectors` module which defines vector
!! algebra routines.
!! 
!! @author Mark Mansfield (formerly of Stevens Institute) and students
!! @author Walid Keyrouz (walid.keyrouz@nist.gov)
!! @date   Thu Sep 19 19:13:59 2013 EDT
!! 
!! @defgroup zeno_vectors Vector Operations
!! 
!! Groups routines related to vector algebra.
!! 
!! @todo
!! Check if Fortran intrinsic functions already perform some of these
!! operations.
!! 
! Time-stamp: <2015-01-05 16:45:33 walid>
! 
! ================================================================

module zeno_vectors

  !! ================================================================

  use numeric_constants

  implicit none
  private

  public normalize, dotproduct, xyzrot, rotate, vector_difference, &
       & cross_product, scalar_product, pythag0, pythag, &
       & backtransform, makepolar

contains

  !! ================================================================
  !! 
  !! Public routines
  !! 
  !! ================================================================

  !> Normalizes vector `n1` and returns result in `n`

  subroutine normalize(n1, n)

    !! Declare arguments

    real, dimension(3) :: n1, n

    !! Local variables

    integer :: i
    real    :: z

    call dotproduct(n1,n1,z)
    z = sqrt(z)
    do i = 1,3
       n(i) = n1(i)/z
    end do

  end subroutine normalize

  !! ----------------------------------------------------------------

  !> @f$ c = a . b @f$, returns value in @c
  !! 
  !! @warning may be replaced by intrinsic routine in Fortran

  subroutine dotproduct(a, b, c)

    !! Declare arguments

    real, dimension(3) :: a,b
    real, intent(out)  :: c

    !! Local variable

    integer :: i

    c = 0.0
    do i = 1,3
       c = c + a(i)*b(i)
    end do

  end subroutine dotproduct

  !! ----------------------------------------------------------------

  !> Generates a 3x3 rotation matrix, `t`, out of three vectors, `n1`,
  !! `n2`, and `n3`.  Returns rotation matrix in `t`
  !! 
  !! @warning does no checking and assumes n1--n3 are unit vectors.

  subroutine xyzrot(n1, n2, n3, t)

    !! Declare arguments

    real, dimension(3)    :: n1, n2, n3
    real, dimension(3, 3) :: t

    !! Local variables

    integer :: i

    do i = 1,3
       t(1,i) = n1(i)
       t(2,i) = n2(i)
       t(3,i) = n3(i)
    end do

  end subroutine xyzrot

  !! ----------------------------------------------------------------

  !> Return in `a` the product a = t.b
  !! 
  !! @param t -- 3x3 real rotation matrix
  !! @param a -- real vector(3)
  !! @param b -- real vector(3)

  subroutine rotate(a, t, b)

    !! Declare arguments

    real, dimension(3)   :: a, b
    real, dimension(3,3) :: t

    !! Local variabls

    integer :: i, j

    do i = 1,3
       a(i) = 0.0
       do j = 1,3
          a(i) = a(i) + t(i,j)*b(j)
       end do
    end do

  end subroutine rotate

  !! ----------------------------------------------------------------

  !> Perform or compute the inverse rotation.  Since it is the reverse
  !! rotation, use the transpose

  subroutine backtransform(c, t, p1)

    !! Declare arguments
    real, dimension(3)   :: c   !> TBD
    real, dimension(3,3) :: t   !> TBD
    real, dimension(3)   :: p1  !> TBD

    !! Local variables
    integer :: i, j
    real, dimension(3) :: v

    do i = 1,3
       v(i) = 0.0
       do j = 1,3
          v(i) = v(i) + t(j,i)*p1(j)
       end do
    end do

    do i = 1,3      
       p1(i) = v(i) + c(i)
    end do

  end subroutine backtransform

  !! ----------------------------------------------------------------

  !> @f$ c = a - b @f$

  subroutine vector_difference(a, b, c)

    !! Declare arguments

    real, dimension(3) :: a, b, c

    !! Local variables

    integer :: i

    do i = 1,3
       c(i) = a(i) - b(i)
    end do

  end subroutine vector_difference

  !! ----------------------------------------------------------------

  !> @f$ c = a \times b @f$

  subroutine cross_product(a, b, c)

    !! Declare arguments

    real, dimension(3) :: a, b, c

    c(3) = a(1)*b(2) - a(2)*b(1)
    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)

  end subroutine cross_product

  !! ----------------------------------------------------------------

  !> c = a_scalar * b_vector; returns result in @c

  subroutine scalar_product(a_scalar, b_vector, c)

    !! Declare arguments

    real, dimension(3) :: b_vector, c
    real, intent(in)   :: a_scalar

    c(1) = a_scalar*b_vector(1)
    c(2) = a_scalar*b_vector(2)
    c(3) = a_scalar*b_vector(3)

  end subroutine scalar_product

  !! ----------------------------------------------------------------

  !> Compute the distance from `x` to the origin by the pythagorean
  !! theorem

  subroutine pythag0(x, d)

    !! Declare arguments

    real, dimension(3) :: x
    real, intent(out)  :: d

    d = sqrt(x(1)**2 + x(2)**2 + x(3)**2)
	
  end subroutine pythag0

  !! ----------------------------------------------------------------

  !> Compute distance between points `x` and `y` by the pythagorean
  !! theorem

  subroutine pythag(x, y, d)

    !! Declare arguments

    real, dimension(3) :: x, y
    real, intent(out)  :: d

    !! Local variables

    integer :: i

    d = 0.0
    do i = 1,3
       d = d + (x(i)-y(i))**2
    end do
    d = sqrt(d)
  
  end subroutine pythag

  !! ----------------------------------------------------------------

  !> Convert unit vector in `p` to spherical-polar coordinates.
  !! 
  !! @warning uses anonymous numeric constants.

  subroutine makepolar(p, theta, phi)

    !! Declare arguments
    real, dimension(3) :: p          !> Unit vector
    real, intent(out)  :: theta, phi !> Spherical coordinates (r == 1)

    if (p(3) < -1.0) then
       theta = M_PI_sp
       phi = 0.0
    else if (p(3) > 1.0) then
       theta = 0.0
       phi = 0.0
    else
       theta = acos(p(3))
       phi = atan2(p(2),p(1))
    end if

  end subroutine makepolar

  !! ================================================================

end module zeno_vectors

!!! ================================================================

!!! Local Variables:
!!! mode: f90
!!! time-stamp-line-limit: 30
!!! fill-column: 72
!!! End:
