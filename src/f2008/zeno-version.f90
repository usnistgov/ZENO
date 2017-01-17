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
!! This file provides the `zeno_version` module which exports constants
!! that describe the version number (Major and Minor) and string.
!! 
!! @author Walid Keyrouz (walid.keyrouz@nist.gov)
!! @date   Tue Feb 24 17:22:46 2015 EST
!! 
!! @defgroup zeno_version Version Number
!! 
!! Defines version number (Major and Minor) and String.
! 
! Time-stamp: <2015-02-25 10:02:39 wtk>
!
! ================================================================

module zeno_version

  !! ================================================================

  implicit none

  public

  integer, parameter :: VERSION_MAJOR = 4
  integer, parameter :: VERSION_MINOR = 0
  character(len=4), parameter :: VERSION_STRING = '4.00'

  !! ================================================================

end module zeno_version

! ================================================================

!!! Local Variables:
!!! mode: f90
!!! time-stamp-line-limit: 30
!!! fill-column: 72
!!! End:
