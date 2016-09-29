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
!! This file defines the `numeric_constants` module defines *named
!! constants* for the rest of the program.
!!
!! @defgroup numeric_constants Numeric Constants
!! 
!! Defines named numeric constants such as `PI` (in single and double
!! precision), etc.
!!
!! @details
!! Module defines well-known single and double precision variants of `PI`
!! (@f$\pi@f$) and @f$e@f$ (the base of natural logarithms).
!!
!! @warning
!! Such constants are not necessarily defined in Fortran 2008.
!!
!! @author Walid Keyrouz (walid.keyrouz@nist.gov)
!! @date   Tue Apr 22 15:00:01 2014 EDT
!
! Time-stamp: <2015-01-21 14:26:56 wtk>
! 
! ================================================================

!> @brief
!! This module defines well-known numeric constants such as `PI`, in both
!! single and double precision floating point numbers, for the rest of
!! the program to avoid repeating literal values in the code.  Their
!! values were copied from `/usr/include/math.h`

module numeric_constants

  !! ================================================================

  use numeric_kinds

  implicit none

  public

  real(dp_k), parameter :: M_PI_dp = 3.14159265358979323846_dp_k
  real(sp_k), parameter :: M_PI_sp = real(M_PI_dp)
  real(dp_k), parameter :: M_PI = M_PI_dp

  real(dp_k), parameter :: M_E_dp = 2.718281828459045235360287471352662498_dp_k
  real(sp_k), parameter :: M_E_sp = real(M_E_dp)
  real(dp_k), parameter :: M_E    = M_E_dp

  !! ================================================================

end module numeric_constants

! ================================================================

!!! Local Variables:
!!! mode: f90
!!! time-stamp-line-limit: 30
!!! fill-column: 72
!!! End:
