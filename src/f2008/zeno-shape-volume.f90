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
!! This file provides the `zeno_shape_volume` module which defines
!! routines that compute the volumes of shapes.
!! 
!! @warning Need definitions for all shapes
!! 
!! @todo Code crying out for an object-oriented aproach with a `shape`
!! class hierarchy.
!!
!! @author Mark Mansfield (formerly of Stevens Institute) and students
!! @author Walid Keyrouz (walid.keyrouz@nist.gov)
!! @date   Fri Oct 25 11:07:40 2013 EDT
!! 
!! @defgroup zeno_shape_volume Shape Volumes
!! 
!! Groups routines that compute shape volumes.
!! 
!! @warning Module is incomplete as it does not handle many shapes.
!! 
!! @warning should be merged as methods in a `shape` class hierarchy.
!! 
! Time-stamp: <2015-01-05 16:45:11 walid>
! 
! ================================================================

module zeno_shape_volume

  !! ================================================================

  use numeric_constants

  implicit none
  private

  public primvol

contains

  !! ================================================================
  !!
  !! Public routines
  !!
  !! ================================================================

  !> Computes the volume of a primitive shape.
  !! 
  !! @warning routine is incomplete!
  !! 
  !! @warning uses anonymous numeric constants.

  subroutine primvol(myelt, pass_vec, tumble, result, early, mess)

    !!GCC$ ATTRIBUTES unused :: tumble

    use zeno_codes_data

    !! Declare arguments
    integer, intent(in)  :: myelt     !> Element type
    real, dimension(12)  :: pass_vec  !> Element specs
    real, dimension(3,3) :: tumble    !> Special rotation matrix
    real, intent(out)    :: result    !> Result returned
    logical, intent(out) :: early     !> TBD
    character(len=10)    :: mess      !> TBD

    !! Local variables
    real, dimension(3) :: x, y, zl, zh
    real :: aa, al, bb, cc, pi, r, r1, r2, side

    !! Get around compiler warning for unused dummy argument(s)
    real, dimension(3,3) :: tumble_2

    !! HACK !!
    if (.false.) then
       tumble_2(1,1) = tumble(1,1)
    end if

    pi = M_PI_sp

    if (myelt == cube_code) then
       side = pass_vec(4)
       result = side**3
       mess = 'CUBE      '
       early = .false.

    else if (myelt == pillar_code) then
       x(1) = pass_vec(1)
       x(2) = pass_vec(5)
       x(3) = pass_vec(9)
       y(1) = pass_vec(2)
       y(2) = pass_vec(6)
       y(3) = pass_vec(10)
       zl(1) = pass_vec(3)
       zl(2) = pass_vec(7)
       zl(3) = pass_vec(11)
       zh(1) = pass_vec(4)
       zh(2) = pass_vec(8)
       zh(3) = pass_vec(12)
       call pillarvol(x,y,zl,zh,result)
       mess = 'PILLAR    '
       early = .false.

    else if (myelt == sphere_code) then
       r = pass_vec(4)
       result = (4.0*pi/3.0)*(r*r*r)
       mess = 'SPHERE    '
       early = .false.

    else if (myelt == triangle_code) then
       early = .true.
       mess = 'TRIANGLE  '

    else if (myelt == disk_code) then
       early = .true.
       mess = 'DISK      '

    else if (myelt == open_cylinder_code) then
       early = .true.
       mess = 'O_CYLINDER'

    else if (myelt == solid_cylinder_code) then
       r = pass_vec(7)
       al = pass_vec(8)
       result = 4.0*pi*r*r*al
       early = .false.
       mess = 'S_CYLINDER'

    else if (myelt == donut_code) then
       r1 = pass_vec(7)
       r2 = pass_vec(8)
       result = 2.0*pi*pi*r1*r2*r2
       early = .false.
       mess = 'TORUS     '

    else if (myelt == ellipsoid_code) then
       aa = pass_vec(10)
       bb = pass_vec(11)
       cc = pass_vec(12)
       result = (4.0*pi/3.0)*aa*bb*cc
       early = .false.
       mess = 'ELLIPSOID '

    else
       stop 'look-up error in primvol'
    end if
        
  end subroutine primvol

  !! ----------------------------------------------------------------

  !> Computes volume of a pillar---Incomplete!
  !! 
  !! @warning Incomplete routine!

  subroutine pillarvol(x, y, zl, zh, result)

    !!GCC$ ATTRIBUTES unused :: x, y, zl, zh, result

    !! Declare arguments
    real, dimension(3) :: x, y, zl, zh    !> Pillar specs
    real, intent(out)  :: result          !> Compute volume returned

    !! Get around compiler warning for unused dummy argument(s)
    real, dimension(3) :: x_2, y_2, zl_2, zh_2

    !! HACK !!
    if (.false.) then
       x_2(1)  = x(1)
       y_2(2)  = y(2)
       zl_2(3) = zl(3)
       zh_2(1) = zh(1)

       result  = 1.0
    end if

    stop 'pillarvol not ready -- sorry'

  end subroutine pillarvol

  !! ================================================================

end module zeno_shape_volume

!!! ================================================================

!!! Local Variables:
!!! mode: f90
!!! time-stamp-line-limit: 30
!!! fill-column: 72
!!! End:
