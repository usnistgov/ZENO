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
!! This file provides the `zeno_rng` module which defines the RNG-related
!! functionality needed in ZENO.
!! 
!! @author Mark Mansfield (formerly of Stevens Institute) and students
!! @author Walid Keyrouz (walid.keyrouz@nist.gov)
!! @date   Tue Sep 10 09:03:13 2013 EDT
!! 
!! @defgroup zeno_rng Random Number Generation
!! 
!! Groups together the three ZENO functions and subroutines that relate
!! to random number generation.  These are the subroutines, `seeder`
!! and `charge`, and the three functions, `ran2`, and `jrand`.
!!
!! Two of these routines, `seeder` and `ran2`, share access to
!! variable `idum`, defined in module `zeno_rng_data`, and store in it
!! the state of the random number generator (RNG).  Several notes:
!!
!! 1. The RNG state is stored in a @b signed @b integer.  Therefore,
!!    the RNG has a period of at most 2^31 = 2,147,483,648 ~= 2E9
!!    which may not be big enough when the number of paths generated
!!    is on the order of 100E6 paths.
!! 
!!    Assuming that each path consists, on average, of 5 points in 3D
!!    space, the set of generated numbers has a cardinality that is
!!    too close to the RNG's period, which may negatively impact the
!!    accuracy of the simulation!
!! 
!! 2. The `jrand` used to declare the `random` common block, but
!!    did not access the variable in it.  Therefore, the declaration
!!    was removed.
!! 
! Time-stamp: <2015-01-25 13:17:44 wtk>
! 
! ================================================================

module zeno_rng

  !! ================================================================

  use fgsl
  use zeno_enums
  use zeno_options, only : rng_opt
  implicit none

  !! ran2 now private
  public seeder, get_rand, jrand, &
       & charge

  private

  integer, save :: idum

  type(fgsl_rng) :: rng_handle
  integer(fgsl_long) :: rng_seed

contains

  !! ================================================================
  !!
  !! Public routines
  !!
  !! ================================================================

  !! Subroutine `seeder` uses a `string` encoding the output of the
  !! `date(1)` command to initialize a seed for the sequence of random
  !! numbers.
  !! 
  !! @param dateline (in) `String` encoding the output of `date(1)`;
  !! assumed to be at least 28 characters long.  An example of a
  !! parameter is @c "Sat Oct 29 14:39:43 MDT 2011".
  !!
  !! @todo Eventually modify so it takes a seed value directly instead
  !! of extracting the seed value from the date-time string.
  !! 
  !! @see nr_ran2()
  !! 
  !! @warning State of RNG is limited to 31-bits.  This limitation is
  !! there from when ZENO used Numerical Recipes functions.

  subroutine seeder(dateline)

    ! This routine uses the UNIX date command to generate a seed for
    ! the sequence of random numbers

    ! Assume 28 chars is minimum size.
    character(len=*), intent(in) :: dateline

    character(len=8) :: shore

    type(fgsl_rng_type) :: t

    if (len(dateline) < 28) then
       stop "Incorrect size of 'dateline'"
    end if

    ! Sat Oct 29 14:39:43 MDT 2011
    ! 0000000001111111111222222222 
    ! 1234567890123456789012345678

    shore(1:2) = dateline(9:10)
    shore(3:4) = dateline(12:13)
    shore(5:6) = dateline(15:16)
    shore(7:8) = dateline(18:19)
    read(shore, *) idum
    idum = -idum

    select case( rng_opt )

    case (RNG_fgsl) 
       t = fgsl_rng_env_setup()
       t = fgsl_rng_default
       rng_handle = fgsl_rng_alloc(t)
       rng_seed = -idum
       call fgsl_rng_set(rng_handle, rng_seed)
    
    case default
       write(*,*) 'Warning: unknown rng option:', rng_opt
       write(*,*) 'Defaults to RNG from (f)gsl'
       t = fgsl_rng_env_setup()
       t = fgsl_rng_default
       rng_handle = fgsl_rng_alloc(t)
       rng_seed = -idum
       call fgsl_rng_set(rng_handle, rng_seed)

    end select

  end subroutine seeder

  !****************************************************************

  !> Simple wrapper to invoke one of the underlying RNGs.

  real function get_rand()

    use zeno_debug, only : incr_RNG_log

#if ! defined(ZENO_ENABLE_LIB_OPTS)

    get_rand = sngl(fgsl_rng_uniform(rng_handle))

#else

    select case( rng_opt )

    case (RNG_fgsl)
       get_rand = sngl(fgsl_rng_uniform(rng_handle))

    case default
       write(*,*) 'Warning: unknown rng option:', rng_opt
       write(*,*) 'Defaults to RNG from (f)gsl'
       get_rand = sngl(fgsl_rng_uniform(rng_handle))

    end select

#endif

    !! Update log data
    call incr_RNG_log()

  end function get_rand

  !> Function `jrand` is a simple wrapper around `ran2` function to
  !! return an integer result.
  !! 
  !! @param k1 (in) integer, lower bound of integer interval for random
  !! integers.
  !! 
  !! @param k2 (in) integer, upper bound of integer interval for random
  !! integers.
  !! 
  !! @return random integer between \p k1 and \p k2, assumes (k1 <= k2).

  integer function jrand(k1, k2)

    integer, intent(in) :: k1, k2

    real :: rr

    if (k2 < k1) stop 'args out of order in jrand'
    if (k1 < 0) stop 'neg. arg in jrand'

    if (k2 == k1) then
       jrand = k1
    else
       do
          !! rr = ran2()
          rr = get_rand()
          rr = float(k2-k1+1)*rr + float(k1)
          jrand = int(rr)
          if (jrand >= k1 .and. jrand <= k2) exit
       end do

    end if

  end function jrand

  !! ----------------------------------------------------------------

  !> Assign three charges to random walker, kk(i) = i-th charge
  !! depending on the RNG and the ratio between rt(i) and r1.
  !! 
  !! @param rt (in)  real vector (3)
  !! @param r1 (in)  real
  !! @param kk (out) integer(3) encoding +1 & -1 random steps
  !! 
  !! @warning rename to something like charges_3_random?

  subroutine charge(rt, r1, kk)

    !! Declare arguments
    real, dimension(3)    :: rt
    real, intent(in)      :: r1
    integer, dimension(3) :: kk

    !! Local variables
    integer :: i
    real    :: p, x

    do i = 1,3
       x = rt(i)/r1
       p = 0.5*(1.0 + x)
       !! if (ran2() < p) then
       if (get_rand() < p) then
          kk(i) = +1
       else
          kk(i) = -1
       end if
    end do

  end subroutine charge

  !! ================================================================

end module zeno_rng

!!! ================================================================

!!! Local Variables:
!!! mode: f90
!!! time-stamp-line-limit: 30
!!! fill-column: 72
!!! End:
