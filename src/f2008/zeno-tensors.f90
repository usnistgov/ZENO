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
!! This file provides the `zeno_tensors` module which defines routines
!! that manipulate tensors.
!! 
!! @author Mark Mansfield (formerly of Stevens Institute) and students
!! @author Walid Keyrouz (walid.keyrouz@nist.gov)
!! @date   Fri Sep 27 16:53:02 2013 EDT
!! 
!! @defgroup zeno_tensors Tensors
!! 
!! Groups routines that manipulate tensors.  It defines the following
!! entry points: `ttdiag`
!! 
! Time-stamp: <2015-02-23 11:10:28 wtk>
!
! ================================================================

module zeno_tensors

  !! ================================================================

  use numeric_kinds
  use zeno_eigen

  implicit none
  private
  public ttdiag, pade

contains

  !! ================================================================
  !!
  !! Public routines
  !!
  !! ================================================================

  !> Diagonalizes the (Rg2 = T) tensor
  !! 
  !! @param tten   -- 3x3 double precision real tensor matrix
  !! @param eigens -- double precision real vector(3) of eigen values

  subroutine ttdiag(tten, eigens)

    !! Declare arguments
    real(dp_k), dimension(3,3) :: tten
    real(dp_k), dimension(3)   :: eigens

    call eigen_vals_3x3_drs_dr(tten, eigens)

  end subroutine ttdiag

  !! ----------------------------------------------------------------

  !> Does the @e Pade approximant to estimate the [eta]/[sigma] ratio,
  !! i.e., q2, from the polarizability tensor.

  subroutine pade(alpha_bongo, q2pade, eigens, xx)

    !! Declare arguments
    real,       dimension(3,3) :: alpha_bongo
    real(dp_k), dimension(3)   :: eigens
    real                       :: q2pade
    real(dp_k), dimension(2)   :: xx

    !! Local variables
    character(len=43), dimension(6,2) :: jet
    integer :: i, j, j1, j2, j3
    integer, dimension(6,2,4) :: var, varmag
    real :: m, alpha, c, bb, b, delta
    real :: num, den
    real :: coef, poww, x1, x2
    real, dimension(6) :: fv

    real, dimension(3) :: d

    call eigen_vals_3x3_srs_dr(alpha_bongo, eigens)

    do i=1,3
       d(i) = real(eigens(i))
    end do

    !> WK: Very bizzare.
    !! 
    !! WK: Why aren't the elts of arrays `var` and `varmag`
    !!     simply initialized instead of read from the `string`
    !!     array, `jet`?
    jet(1,1) = 'delta 4800 -3    660 -3  -1247 -3    787 -3'
    jet(1,2) = 'k        0 -3   1040 -3   2012 -3   2315 -3'
    jet(2,1) = 'b      680 -3  -7399 -3   1048 -3    136 -3'
    jet(2,2) = 't        0 -3   1063 -3    895 -3   4993 -3'
    jet(3,1) = 'B     1925 -3  -8611 -3   1652 -3   -120 -3'
    jet(3,2) = 'q        0 -3   1344 -3   2029 -3   1075 -3'
    jet(4,1) = 'c     1343 -2   1617 -2     51 -2   -586 -2'
    jet(4,2) = 'r        0 -3    489 -3    879 -3   2447 -3'
    jet(5,1) = 'alpha 1623 -2  -1592 -2   1483 -2   -374 -2'
    jet(5,2) = 'v        0 -3    462 -3   1989 -3   4608 -3'
    jet(6,1) = 'm     2786 -3    293 -3   -110 -3     12 -3'
    jet(6,2) = 'u        0 -3    556 -3   2034 -3   3024 -3'
    !>          1234567890123456789012345678901234567890123

    do i = 1,6
       do j = 1,2
          jet(i,j)(1:5) = '     '
       end do
    end do

    do j1 = 1,6
       do j2 = 1,2
          read(jet(j1,j2),*) (var(j1,j2,j3),varmag(j1,j2,j3),j3=1,4)
       end do
    end do

    !> @warning compares two reals with ==
    if (d(1) /= 0.0) then

       xx(1) = dble(d(2)/d(1))
       xx(2) = dble(d(3)/d(2))
       x1 = alog(d(2)/d(1))
       x2 = alog(d(3)/d(2))

       do j1 = 1,6
          fv(j1) = 0.0
          do j3 = 1,4
             coef = float(var(j1,1,j3)) * 10.0**(varmag(j1,1,j3))
             poww = float(var(j1,2,j3)) * 10.0**(varmag(j1,2,j3))
             fv(j1) = fv(j1) + coef*exp(-poww*x1)
          end do
       end do

    else

       !> The object must be two-dimensional
       !> x1 = infinity

       xx(1) = dble(d(2)/d(1))
       xx(2) = dble(d(3)/d(2))
       x2 = alog(d(3)/d(2))
       do j1 = 1,6
          coef = float(var(j1,1,1)) * 10.0**(varmag(j1,1,1))
          fv(j1) = coef
       end do

    end if

    delta = fv(1)
    b = fv(2)
    bb = fv(3)
    c = fv(4)
    alpha = fv(5)
    m = fv(6)

    num = delta*alpha + c*x2 + b*x2**2 + 4.0*(x2**m)
    den = 6.0*alpha + (6.0*c*x2/delta) + bb*x2**2 + 5.0*(x2**m)
    q2pade = num/den

  end subroutine pade

  !! ================================================================

end module zeno_tensors

!!! ================================================================

!!! Local Variables:
!!! mode: f90
!!! time-stamp-line-limit: 30
!!! fill-column: 72
!!! End:
