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
!! This file provides the `numeric_kinds` module which defines numeric 
!! **kinds** for the rest of the program.
!! 
!! @defgroup numeric_kinds Numeric Kinds
!! 
!! Defines **numeric kinds**, with different sizes and precisions, for
!! the rest of the program.
!! 
!! @details
!! Module defines numeric kinds (`integer` and `real`) for the rest of
!! the program.  It defines:
!! - 3 signed integers of different sizes (2, 4, and 9) and
!! - 3 reals (sp, dp, and qp) with different precisions.
!! 
!! Sources:
!! @verbatim 
!!   p. 71
!!   Metcalf, Reid & Cohen
!!   Modern Fortran explained
!!   Oxford University Press, (c) 2011
!! @endverbatim
!! 
!! @see
!! - http://fortranwiki.org/fortran/show/Real+precision
!! - http://fortranwiki.org/fortran/show/selected_int_kind
!! 
!! @author Walid Keyrouz (walid.keyrouz@nist.gov)
!! @date   Thu Sep 19 14:37:07 2013 EDT
! 
! Time-stamp: <2015-01-20 18:40:44 wtk>
! 
! ================================================================

module numeric_kinds

  !! ================================================================

  use, intrinsic :: iso_fortran_env

  implicit none

  !> @note
  !! Named constants, `iSb_k`, stand for `S = 1, 2, 4, 8` byte
  !! integers param @f$N@f$ of `selected_int_kind` specifies integer
  !! range: @f$] -10^N, +10^N [@f$
  !! 
  !! @see http://fortranwiki.org/fortran/show/selected_int_kind

  integer, parameter :: i1b_k = selected_int_kind(2)  !< 1-byte signed int
  integer, parameter :: i2b_k = selected_int_kind(4)  !< 2-byte signed int
  integer, parameter :: i4b_k = selected_int_kind(9)  !< 4-byte signed int
  integer, parameter :: i8b_k = selected_int_kind(18) !< 8-byte signed int

  ! and for single, double, and quadruple precision reals:
  ! ref: http://fortranwiki.org/fortran/show/Real+precision

  integer, parameter :: sp_k = REAL32     !< Single precision
  integer, parameter :: dp_k = REAL64     !< Double precision
  integer, parameter :: qp_k = REAL128    !< Quad precision

  ! Old defs from Metcalf, Reid, & Cohen's book
  ! integer, parameter :: sp_k = kind(1.0)
  ! integer, parameter :: dp_k = selected_real_kind(2*precision(1.0_sp_k))
  ! integer, parameter :: qp_k = selected_real_kind(2*precision(1.0_dp_k))

end module numeric_kinds

!!! ================================================================

!!! Local Variables:
!!! mode: f90
!!! time-stamp-line-limit: 30
!!! fill-column: 72
!!! End:
