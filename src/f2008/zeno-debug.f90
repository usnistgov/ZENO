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
!! This file defines the `zeno_debug` module which groups together
!! routines that generate `debug` and `log` output.
!! 
!! @author Walid Keyrouz (walid.keyrouz@nist.gov)
!! @date   Tue Dec 24 15:27:52 2013 EST
!! 
!! @defgroup zeno_debug Debugging and Logging Utilities
!! 
!! Groups functions & routines that generate `debug` and `log` output.
!! There are two compile-time flags that may be enabled:
!! 
!! - ZENO_ENABLE_RNG_LOG -- enables counting the RNG calls.

!! - ZENO_ENABLE_LIB_OPTS -- enables options to specify libraries for
!!   (1) generating random numbers (currently only GSL is available) and
!!   (2) performing eigen-value computations.
!! 
!! @warning Possibly use a more general `logging` and 
!! `error-reporting` mechanism.
!
! Time-stamp: <2015-01-05 16:43:07 walid>
! 
! ================================================================

module zeno_debug

  !! ================================================================

  use numeric_kinds
  use zeno_enums

  implicit none

  public write_RNG_log,  print_RNG_log, incr_RNG_log, &
       & write_LS_log,   print_LS_log, &
       & write_DIST_log, print_DIST_log

  !> Controls debug output stmnts.
  logical, public, save :: silent = .false.

  !> Controls `curses` tricks.
  logical, public, save :: print_flush = .false.

  !> Calls to `listersort`
  integer(i8b_k), public, save :: calls_LS_cnt = 0

  !> Calls to @ distance
  integer(i8b_k), public, save :: calls_DIST_cnt = 0

  !> Calls to `distance_scan_neighbors`
  integer(i8b_k), public, save :: calls_DIST_scan_neighbors = 0

  !> Calls to `distance_scan_all_elts`
  integer(i8b_k), public, save :: calls_DIST_scan_all_elts = 0

  !> Calls to `re_compute_bubble`
  integer(i8b_k), public, save :: calls_DIST_re_comp_bbl = 0

  private

  !> Calls to RNG
  integer(i8b_k), save :: calls_rng = 0

contains

  !! ================================================================
  !!
  !! Public routines
  !!
  !! ================================================================

  !> Prints the number of times a random number is generated.
  !! 
  !! @param funit  I/O unit to be used for @e logging output

  subroutine write_RNG_log(funit)

    integer, intent(in) :: funit


#if defined(ZENO_ENABLE_RNG_LOG)

    write(funit, '(''RNG #calls: '', I0)') calls_rng

#else

    write(funit, '(a)') 'ZENO must be compiled with ZENO_ENABLE_RNG_LOG'

#endif

  end subroutine write_RNG_log

  !! ----------------------------------------------------------------

  !> Prints the number of times a random number is generated to
  !! Fortran's equivalent of `stdout`.

  subroutine print_RNG_log()

    use, intrinsic :: iso_fortran_env, only : OUTPUT_UNIT

    call write_RNG_log(OUTPUT_UNIT)

  end subroutine print_RNG_log

  !! ----------------------------------------------------------------

  !> Increments the number of calls to generate a random number

  subroutine incr_RNG_log()

#if defined(ZENO_ENABLE_RNG_LOG)

    calls_rng = calls_rng + 1

#endif

  end subroutine incr_RNG_log

  !! ----------------------------------------------------------------

  !> Prints the call count of `listersort`.
  !! 
  !! @param funit  I/O unit to be used for @e logging output

  subroutine write_LS_log(funit)

    !! Declare arguments
    integer :: funit

    write (funit, '(''LS #calls: '', I0)') calls_LS_cnt

  end subroutine write_LS_log

  !! ----------------------------------------------------------------

  !> Prints the call count of `listersort` to Fortran's equivalent of
  !! `stdout`.

  subroutine print_LS_log()

    use, intrinsic :: iso_fortran_env, only : OUTPUT_UNIT

    call write_LS_log(OUTPUT_UNIT)

  end subroutine print_LS_log

  !! ----------------------------------------------------------------

  !> Prints the call counts of `distance`, `distance_scan_neighbors`,
  !! and `distance_scan_all_elts`.
  !! 
  !! @param funit  I/O unit to be used for @e logging output

  subroutine write_DIST_log(funit)

    !! Declare arguments
    integer :: funit

    write(funit, '(''DIST #calls: '', I0)') calls_DIST_cnt
    write(funit, '(''DIST scan neighbors (bbl_rad >= 0): '', I0)') &
         & calls_DIST_scan_neighbors
    write(funit, '(''DIST scan all elts: '', I0)') calls_DIST_scan_all_elts
    write(funit, '(''DIST re comp bubble (pd+ds > bbl_rad): '', I0)') &
         & calls_DIST_re_comp_bbl

  end subroutine write_DIST_log

  !! ----------------------------------------------------------------

  !> Prints the call counts of `distance`, `distance_scan_neighbors`,
  !! and `distance_scan_all_elts` to Fortran's equivalent of `stdout`.

  subroutine print_DIST_log()

    use, intrinsic :: iso_fortran_env, only : OUTPUT_UNIT

    call write_DIST_log(OUTPUT_UNIT)

  end subroutine print_DIST_log

  !! ================================================================

end module zeno_debug

!!! ================================================================

!!! Local Variables:
!!! mode: f90
!!! time-stamp-line-limit: 30
!!! fill-column: 72
!!! End:
