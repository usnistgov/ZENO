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
!! This file provides the zeno_options module which defines some global
!! options used throughout ZENO.  This was put place while evolving ZENO
!! to remove Numerical Recipes code from the code base.
!! 
!! @author Walid Keyrouz (walid.keyrouz@nist.gov)
!! @date   Tue May  6 12:10:56 2014 EDT
!! 
!! @defgroup zeno_options Global Options
!! 
!! Defines some global options used throughout ZENO.
!! 
!! @warning module may disappear altogether.
!! 
! Time-stamp: <2015-01-26 12:58:30 wtk>
! 
! ================================================================

module zeno_options

  !! ================================================================

  use zeno_enums

  implicit none

  !> The default choice for eigenvalue computation routine
  integer, public, save :: eig_opt = EIG_wp

  !> The default choice for the random number generator
  integer, public, save :: rng_opt = RNG_fgsl

  !! ================================================================

end module zeno_options

!!! ================================================================

!!! Local Variables:
!!! mode: f90
!!! time-stamp-line-limit: 30
!!! fill-column: 72
!!! End:
