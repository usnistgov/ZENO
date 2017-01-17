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
!! This file provides the `zeno_sphere` module which defines routines
!! that relate to spheres.
!! 
!! @author Mark Mansfield (formerly of Stevens Institute) and students
!! @author Walid Keyrouz (walid.keyrouz@nist.gov)
!! @date   Fri Sep 20 15:01:19 2013 EDT
!! 
!! @defgroup zeno_sphere Sphere
!! 
!! Exports sphere routine only.  Defined in a separate module because
!! the routine uses random number generators and is used by the @c
!! zeerot routine!
!!
!! @warning Liable to be changed.
!!
!! @warning Functionality should be generalized and integrated into a @c
!! shape class hierarchy.
! 
! Time-stamp: <2015-01-26 13:15:52 wtk>
! 
! ================================================================

module zeno_sphere

  !! ================================================================

  use zeno_rng, only : get_rand

  implicit none

  public sphere

contains

  !! ================================================================
  !!
  !! Public routine
  !!
  !! ================================================================

  !> Generate a point distributed randomly on a spherical surface of
  !! radius `r1`
  !! 
  !! @param rt: real array of size 3 to be filled
  !! @param r1: radius of sphere
  !! 
  !! @warning Consider changing algorithm to one that uses fewer random
  !! numbers.
  !! 
  !! @todo Rename routine to something like `point_on_sphere` or 
  !! `point_on_shape`.

  subroutine sphere(rt, r1)

    !! Declare arguments

    real, dimension(3) :: rt    !> unit vector generated
    real, intent(in)   :: r1    !> sphere radius

    !! Local variables

    integer :: j
    real    :: dd

    do
       do j = 1,3
          rt(j) = get_rand()*2.0 - 1.0
       end do

       dd = rt(1)**2 + rt(2)**2 + rt(3)**2

       if (dd <= 1.0) exit
    end do

    dd = sqrt(dd)
    do j = 1,3
       rt(j) = r1*rt(j)/dd
    end do

  end subroutine sphere

end module zeno_sphere

!!! ================================================================

!!! Local Variables:
!!! mode: f90
!!! time-stamp-line-limit: 30
!!! fill-column: 72
!!! End:
