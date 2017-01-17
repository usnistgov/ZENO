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
!! This file provides the `zeno_files_data` module which defines global
!! varialbles that relate to file names and I/O units.
!! 
!! @author Walid Keyrouz (walid.keyrouz@nist.gov)
!! @date   Mon Dec 23 18:29:42 2013 EST
!! 
!! @defgroup zeno_files_data File Names and Units for ZENO files
!! 
!! Replaces `/filenames/` and `/filenumbers/` common blocks.
!!
!! @todo document all data members.
! 
! Time-stamp: <2015-02-24 16:07:19 walid>
! 
! ================================================================

module zeno_files_data

  !! ================================================================

  implicit none

  public

  !> I/O unit for `".bod"` file
  integer, save :: nbod

  integer, save :: &
       &  nzno,nznr,nstk,ndfl,nefl,nzsd,nzh,nih,nsh,nph,nlog

  !> Name for `".bod"` file
  character(len=25), save :: fbod

  character(len=25), save :: &
       &  fzno,fznr,fstk,fdfl,fefl,fzsd,fzh,fih,fsh,fph,flog

  !! ================================================================

end module zeno_files_data

!!! ================================================================

!!! Local Variables:
!!! mode: f90
!!! time-stamp-line-limit: 30
!!! fill-column: 72
!!! End:
