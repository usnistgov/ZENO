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
!! This file provides the `zeno_zeerot` module which generates rotation
!! matrices.
!! 
!! @author Mark Mansfield (formerly of Stevens Institute) and students
!! @author Walid Keyrouz (walid.keyrouz@nist.gov)
!! @date   Fri Sep 20 15:57:19 2013 EDT
!! 
!! @defgroup zeno_zeerot Z-axis Rotation Matrices?
!!
!! Module exports zeerot routine only.  Used because routine uses
!! random number generator and sphere routine to generate rotation
!! matrix.  There is much more direct approach to do so!
!! 
!! @warning Liable to be changed to a far more direct approach such as
!! the one documented by Glenn Murray in
!! http://inside.mides.edu/~gmurray/ArbitraryAxisRotation
!! 
! Time-stamp: <2015-01-05 16:46:07 walid>
! 
! ================================================================

module zeno_zeerot

  !! ================================================================

  use zeno_rng, only : get_rand
  use zeno_vectors, only : dotproduct, cross_product
  use zeno_sphere

  implicit none

  public zeerot

  private
  
  real, parameter :: DOT_PROD_MIN = 0.1

contains

  !! ================================================================
  !!
  !! Public routine
  !!
  !! ================================================================

  !> Returns a rotation matrix in `t` that rotates `n3` into the
  !! positive z-axis.
  !! 
  !! @param n3 -- real unit vector(3)
  !! @param t -- filled 3x3 real rotation matrix
  !! 
  !! @pre `n3` is a unit vector; does no checking
  !! @post `t` is a 3x3 rotation matrix; trans(t) = inv(t), etc.
  !! 
  !! @warning
  !! Odd way of generating rotation matrix.  Here it is: Proceed by
  !! taking an arbitrary normalized vector, but throw it away and start
  !! over if it is too close to `n3`.  The second axis is taken as a
  !! linear combination of these two vectors to be a unit vector.  The
  !! third one is a cross product of the two.
  !! 
  !! @todo
  !! Consider using a direct approach as documented by Glenn Murray in
  !! http://inside.mines.edu/~gmurray/ArbitraryAxisRotation

  subroutine zeerot(n3, t)

    !! Declare arguments
    real, dimension(3)   :: n3
    real, dimension(3,3) :: t

    !! Local variables

    real, dimension(3) :: n1, n2, nx
    integer            :: i
    real               :: alpha, beta, dx3

    do
       call sphere(nx,1.0)
       call dotproduct(nx,n3,dx3)

       if (1.0-abs(dx3) > DOT_PROD_MIN) exit
    end do

    beta = sqrt(1.0/(1.0 - dx3**2))
    alpha = -beta*dx3
    do i = 1,3
       n1(i) = alpha*n3(i) + beta*nx(i)
    end do
    call cross_product(n3,n1,n2)
 
    do i = 1,3
       t(1,i) = n1(i)
       t(2,i) = n2(i)
       t(3,i) = n3(i)
    end do

  end subroutine zeerot

  !! ================================================================

end module zeno_zeerot

!!! ================================================================

!!! Local Variables:
!!! mode: f90
!!! time-stamp-line-limit: 30
!!! fill-column: 72
!!! End:
