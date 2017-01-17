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
!! This file provides the `zeno_enums` module which defines `enum` types
!! used in ZENO in general and, more particularly, in handling
!! command-line options.
!!
!! @author Walid Keyrouz (walid.keyrouz@nist.gov)
!! @date   Tue Apr 29 16:39:56 2014 EDT
!! 
!! @defgroup zeno_enums Enumerated Types
!! 
!! Defines `enum` types that are used in ZENO in general, and more
!! particularly, in handling command-line options.  For now, the enums
!! are used to define the values of `options`.
!!
!! This module was used to evolve the ZENO code-base by adding code
!! that replaced routines based on Numerical Recipes, 2/e, for
!! eigenvalue computations and random number generation.
!!
! Time-stamp: <2015-01-05 16:43:47 walid>
! 
! ================================================================

module zeno_enums

  !! ================================================================

  use, intrinsic :: iso_c_binding

  implicit none

  !! Choose algorithm for computing eigen values.
  enum, bind(c)
     enumerator :: EIG_wp    = 1, EIG_WikiPedia = 1
  end enum

  !! Choose library to generate random numbers.
  enum, bind(c)
     enumerator :: RNG_fgsl  = 1
  end enum

  !! ================================================================

end module zeno_enums

!!! ================================================================

!!! Local Variables:
!!! mode: f90
!!! time-stamp-line-limit: 30
!!! fill-column: 72
!!! End:
