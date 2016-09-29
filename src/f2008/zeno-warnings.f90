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
!! This file provide the `zeno_warnings` module which replaces the old 
!! `/cubit/` and `/sell/` common blocks.
!! 
!! @author Mark Mansfield (formerly of Stevens Institute) and students
!! @author Walid Keyrouz (walid.keyrouz@nist.gov)
!! @date   Fri Jan 17 16:39:45 2014 EST
!! 
!! @defgroup zeno_warnings Warning Messages
!! 
!! Module declares data members that used to be in `/cubit/` and 
!! `/sell/` common blocks.
!! 
!! - Members `ncube` and `ferr` warn about points that are @e
!!   slightly inside a cube.  The user will be told how far inside
!!   the cube the offending point(s) are found.
!! 
!! - Members `nell` and `rerr` warn about @e overstretching of
!!   ellipsoids.  Coder of ellipsoid overstretching code was not
!!   completely confident of the stretching equations.  If these
!!   members get set, it is an indication that over-stretching has
!!   happened.  The user will be warned about this with an error
!!   statement.
!! 
!! - Member `ncube`, an `integer`
!! 
!! - Member `ferr`, a `real`.
!! 
!! - Member `nell`, an `integer`
!! 
!! - Member `rerr`, a `real`.
!! 
! Time-stamp: <2015-01-05 16:45:38 walid>
! 
! ================================================================

module zeno_warnings

  !! ================================================================

  implicit none

  public
        
  integer :: ncube              !> Number of offending cubes
  real    :: ferr               !> TBD
  integer :: nell               !> Number of offending ellipsoids
  real    :: rerr               !> TBD

  !! ================================================================

end module zeno_warnings

!!! ================================================================

!!! Local Variables:
!!! mode: f90
!!! time-stamp-line-limit: 30
!!! fill-column: 72
!!! End:
