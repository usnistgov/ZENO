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
!! This file defines the `zeno_eigen` module which provides the routines
!! for computing eigen values in ZENO.
!!
!! @author Walid Keyrouz (walid.keyrouz@nist.gov)
!! @date   Fri Apr 18 17:44:58 2014 EDT
!!
!! @defgroup zeno_eigen Eigenvalue Computation Routines
!! 
!! Computes eigenvalues of a 3x3 real symmetric matrix, single or double
!! precision.  Eigenvalues are returned as a sorted vector of 3 doubles.
!!
!! @details
!! Internally uses code listed on Wikipedia.  This was done to minimize
!! dependencies on external libraries.
!! 
!! @see http://en.wikipedia.org/wiki/Eigenvalue_algorithm#3.C3.973_matrices
!
! Time-stamp: <2015-01-23 17:16:16 wtk>
!
! ================================================================

module zeno_eigen

  !! ================================================================

  use numeric_kinds
  use numeric_constants

  use zeno_enums
  use zeno_options

  implicit none

  private
  public eigen_vals_3x3_drs_dr, eigen_vals_3x3_srs_dr

contains

  !! ================================================================
  !!
  !! Public routines
  !!
  !! ================================================================

  !> Computes the eigen values of a 3x3 @b Double precision Real
  !! Symmetric matrix (3x3_drs)
  !! 
  !! @param a (in) 3x3 double precision real array (`3x3_drs`)
  !! @param lambda (out) 3 double precision real vector (`dr`).

  subroutine eigen_vals_3x3_drs_dr(a, lambda)

    !! Declare arguments
    real(dp_k), dimension(3,3) :: a
    real(dp_k), dimension(3)   :: lambda

#if ! defined(ZENO_ENABLE_LIB_OPTS)

    !> Wikipedia 3x3 algorithm
    call WP_eigen_vals_3x3_drs_dr(a, lambda)

#else

    select case( eig_opt )

       case( Eig_WP )              !> Wikipedia 3x3 algorithm
          call WP_eigen_vals_3x3_drs_dr(a, lambda)

       case default
          write(*,*) 'Warning: unknown eigen option:', eig_opt
          write(*,*) 'Defaults to algorithm from WikiPedia, 3x3 cases'
          call WP_eigen_vals_3x3_drs_dr(a, lambda)

    end select

#endif

  end subroutine eigen_vals_3x3_drs_dr

  !! ----------------------------------------------------------------

  !> Computes the eigen values of a 3x3 @b Single precision Real
  !! Symmetric matrix (3x3_drs)
  !! 
  !! @param s_a    (in)  3x3 single precision real array (`3x3_drs`)
  !! @param lambda (out) 3 double precision real vector (`dr`).

  subroutine eigen_vals_3x3_srs_dr(s_a, lambda)

    !! Declare arguments
    real(sp_k), dimension(3,3) :: s_a       !> input matrix
    real(dp_k), dimension(3)   :: lambda    !> output vector

    !! Local vars
    real(dp_k), dimension(3,3) :: d_a

    integer :: i, j

    select case( eig_opt )

       case( Eig_WP )              !> Wikipedia 3x3 algorithm

          do i = 1, 3
             do j =  1, 3
                d_a(i, j) = s_a(i, j)
             end do
          end do

          call WP_eigen_vals_3x3_drs_dr(d_a, lambda)

       case default
          write(*,*) 'Warning: unknown eigen option:', eig_opt
          write(*,*) 'Defaults to WikiPedia, 3x3 cases'

          do i = 1, 3
             do j =  1, 3
                d_a(i, j) = s_a(i, j)
             end do
          end do

          call WP_eigen_vals_3x3_drs_dr(d_a, lambda)

    end select

  end subroutine eigen_vals_3x3_srs_dr


  !! ================================================================
  !! 
  !! Private routines
  !! 
  !! ================================================================

  !> Computes the eigen values of a 3x3 @b Double precision Real
  !! Symmetric matrix (3x3_drs).  The algorithm implemented is the one
  !! listed in Wikipedia at the page given below.
  !! 
  !! @param a (in) 3x3 single precision real array (`3x3_drs`)
  !! @param lambda (out) 3 double precision real vector (`dr`).
  !!
  !! @see http://en.wikipedia.org/wiki/Eigenvalue_algorithm#3.C3.973_matrices

  subroutine WP_eigen_vals_3x3_drs_dr(a, lambda)

    !! Declare arguments
    real(dp_k), dimension(3,3) :: a          !> input matrix
    real(dp_k), dimension(3)   :: lambda     !> output vector

    !! Local variables
    real(dp_k) :: p, p1, p2
    real(dp_k) :: q
    real(dp_k) :: det_b, r, phi
    real(dp_k) :: eig1, eig2, eig3
    real(dp_k), dimension(3,3) :: b

    integer :: i, j             ! array indices
    
    !! Given a real symmetric 3x3 matrix A, compute the eigenvalues
 
    p1 = a(1,2)*a(1,2) + a(1,3)*a(1,3) + a(2,3)*a(2,3)

    if (p1 == 0) then
       !! A is diagonal.
       lambda(1) = a(1,1)
       lambda(2) = a(2,2)
       lambda(3) = a(3,3)

       !! Sort vector
       if (lambda(1) > lambda(2)) then
          eig1 = lambda(1)
          lambda(1) = lambda(2)
          lambda(2) = eig1
       end if

       if (lambda(2) > lambda(3)) then
          eig2 = lambda(2)
          lambda(2) = lambda(3)
          lambda(3) = eig2
       end if

       if (lambda(1) > lambda(2)) then
          eig1 = lambda(1)
          lambda(1) = lambda(2)
          lambda(2) = eig1
       end if

    else

       q = (a(1,1) + a(2,2) + a(3,3)) / 3
       p2 =   (a(1,1)-q)*(a(1,1)-q) + (a(2,2)-q)*(a(2,2)-q) + &
            & (a(3,3)-q)*(a(3,3)-q) + 2*p1
       p = sqrt(p2 / 6)

       b(1,1) = a(1,1) - q
       b(1,2) = a(1,2)
       b(1,3) = a(1,3)

       b(2,1) = a(2,1)
       b(2,2) = a(2,2) - q
       b(2,3) = a(2,3)

       b(3,1) = a(3,1)
       b(3,2) = a(3,2)
       b(3,3) = a(3,3) - q

       do i = 1,3
          do j = 1,3
             b(i,j) = b(i,j) / p
          end do
       end do

       !! Use Sarrus's rule (http://wikipedia.org/wiki/Rule_of_Sarrus)
       det_b =  b(1,1) * b(2,2) * b(3,3) + b(1,2) * b(2,3) * b(3,1) &
            & + b(1,3) * b(2,1) * b(3,2) - b(3,1) * b(2,2) * b(1,3) &
            & - b(3,2) * b(2,3) * b(1,1) - b(3,3) * b(2,1) * b(1,2)

       r = det_b / 2
 
       !! In exact arithmetic for a symmetric matrix  -1 <= r <= 1
       !! but computation error can leave it slightly outside this range.
       if (r <= -1) then
          phi = M_PI / 3
       else if (r >= 1) then
          phi = 0
       else
          phi = acos(r) / 3
       end if
 
       !! the eigenvalues satisfy eig3 <= eig2 <= eig1
       eig1 = q + 2 * p * cos(phi)
       eig3 = q + 2 * p * cos(phi + (2*M_Pi/3))
       eig2 = 3 * q - eig1 - eig3 !! since trace(A) = eig1 + eig2 + eig3

       lambda(1) = eig3
       lambda(2) = eig2
       lambda(3) = eig1
    end if

  end subroutine WP_eigen_vals_3x3_drs_dr

  !! ================================================================

end module zeno_eigen

!!! ================================================================

!!! Local Variables:
!!! mode: f90
!!! time-stamp-line-limit: 30
!!! fill-column: 72
!!! End:
