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
!! This file defines the `zeno_stat_utils` module which groups some
!! statistical utility routines.
!! 
!! @author Mark Mansfield (formerly of Stevens Institute) and students
!! @author Walid Keyrouz (walid.keyrouz@nist.gov)
!! @date   Fri Nov  1 15:47:42 2013 EDT
!!
!! @defgroup zeno_stat_utils Stat Utilities
!! 
!! Groups statistics utility routines used in ZENO.
!! 
!! @warning Routines to be renamed to be more descriptive of
!! functionality.
!! 
!! @todo check if these can be replaced by Fortran intrinsics.
!! 
! Time-stamp: <2015-01-05 16:42:53 walid>
! 
! ================================================================

module zeno_stat_utils

  !! ================================================================

  use numeric_kinds

  implicit none
  private
  public mean20, tally, accume

contains

  ! ================================================================
  !
  ! Public routines
  !
  ! ================================================================

  !> @brief Computes the average and standard deviation (@f$\sigma@f$)
  !! of 20 real values.  `delta` guaranteed to be >= 0.
  !! 
  !! @param xar   (in)  array of @b 20 reals
  !! 
  !! @param x     (out) average
  !! @param delta (out) standard deviation, sigma
  !!
  !! @warning uses a fixed-size array!

  subroutine mean20(xar, x, delta)

    !! Declare arguments
    real, dimension(20) :: xar    !< Array of 20 reals
    real, intent(out)   :: x      !> average
    real, intent(out)   :: delta  !> standard deviation, @f$\sigma@f$

    !! Local variables
    integer :: i
    real    :: s1,s2

    s1 = 0.0
    s2 = 0.0

    do i = 1,20
       s1 = s1 + xar(i)
       s2 = s2 + xar(i)**2
    end do

    x = s1/20.0
    s1 = s1/20.0
    s2 = s2/20.0
    s2 = s2 - s1**2
    if (s2 < 0.0) s2 = 0.0
    delta = sqrt(s2/20.0)

  end subroutine mean20

  ! ----------------------------------------------------------------

  !> Convert raw statistics on random walk trajectories into
  !! polarizability tensor and capacitance
  !!
  !! `tae`, `uae`,`vae`, `wae` are being calculated and reported
  !! for the cases in which an ensemble average will be taken
  !! 
  !! @warning argument list is too long!  

  subroutine tally(khitp, khite, vp, ve, sum, aa, daa, cap, delta_cap, &
       &           r1, tae, uae, vae, wae, rlaunch)

    !!GCC$ ATTRIBUTES unused :: rlaunch

    !! Declare arguments
    integer, dimension(3,20) :: khitp        !> TBD
    integer, dimension(3,20) :: khite        !> TBD
    real(dp_k), dimension(3,3,20) :: vp      !> TBD
    real(dp_k), dimension(3,3,20) :: ve      !> TBD
    real(dp_k), dimension(20)     :: sum     !> TBD
    real, dimension(3,3) :: aa               !> TBD
    real                 :: cap              !> TBD
    real, dimension(3,3) :: daa              !> TBD
    real                 :: delta_cap        !> TBD
    real                 :: r1               !> TBD
    real(dp_k), dimension(3)   :: tae        !> TBD
    real(dp_k), dimension(3)   :: uae        !> TBD
    real(dp_k), dimension(3,3) :: vae        !> TBD
    real(dp_k), dimension(3,3) :: wae        !> TBD
    real, intent(in)           :: rlaunch    !> TBD

    !! Local variables
    real(dp_k) tt(3,20),tu(3,20)
    real(dp_k), dimension(3,3,20) :: tv, tw
    real, dimension(20) :: capar
    real, dimension(3,3,20) :: aar
    real(dp_k) :: nan
    integer :: i, j, k

    !! Get around compiler warning for unused arguments
    real :: rlaunch_2

    !! HACK !!
    if (.false.) then
       rlaunch_2 = rlaunch
    end if

    nan = 0.0d0
    do i = 1,20
       nan = nan + sum(i)
    end do

    do i = 1,3
       do j = 1,20
          tt(i,j) = dble(khitp(i,j)+khite(i,j))/sum(j)
          tu(i,j) = dble(khitp(i,j)-khite(i,j))/sum(j)
       end do
    end do

    do i = 1,3
       tae(i) = 0.0d0
       uae(i) = 0.0d0
       do j = 1,20
          tae(i) = tae(i) + dble(khitp(i,j)+khite(i,j))
          uae(i) = uae(i) + dble(khitp(i,j)-khite(i,j))
       end do
       tae(i) = tae(i)/nan
       uae(i) = uae(i)/nan
    end do

    do i = 1,3
       do j = 1,3
          do k = 1,20
             tv(i,j,k) = (vp(i,j,k) + ve(i,j,k))/sum(k)
             tw(i,j,k) = (vp(i,j,k) - ve(i,j,k))/sum(k)
          end do
       end do
    end do

    do i = 1,3
       do j = 1,3
          vae(i,j) = 0.0d0
          wae(i,j) = 0.0d0
          do k = 1,20
             vae(i,j) = vp(i,j,k) + ve(i,j,k) + vae(i,j)
             wae(i,j) = vp(i,j,k) - ve(i,j,k) + wae(i,j)
          end do
          vae(i,j) = vae(i,j)/nan
          wae(i,j) = wae(i,j)/nan
       end do
    end do

    do k = 1,20
       capar(k) = REAL(tt(1,k)*r1)
    end do

    do i = 1,3
       do j = 1,3
          do k = 1,20
             aar(i,j,k) = sngl(tw(i,j,k)-tu(j,k)*tv(i,j,k)/tt(j,k))
          end do
       end do
    end do

    do i = 1,3
       do j = 1,3
          do k = 1,20
             aar(i,j,k) = 3.0 * r1 * r1 * aar(i,j,k)
          end do
       end do
    end do

    call mean20(capar,cap,delta_cap)

    do i = 1,3
       do j = 1,3
          do k = 1,20
             capar(k) = aar(i,j,k)
          end do
          call mean20(capar,aa(i,j),daa(i,j))
       end do
    end do

  end subroutine tally

  !! ----------------------------------------------------------------

  !> Accumulate statistics on random walks
  !!
  !! It's probably not necessary at present, but just in case this
  !! ever gets run with very large numbers of trajectories, accumulate
  !! statistics in double precision

  subroutine accume(t, kk, khitp, khite, vp, ve, loop)

    !! Declare arguments
    real, dimension(3)            :: t      !> TBD
    integer, dimension(3)         :: kk     !> TBD
    integer, dimension(3,20)      :: khitp  !> TBD
    integer, dimension(3,20)      :: khite  !> TBD
    real(dp_k), dimension(3,3,20) :: vp, ve !> TBD
    integer                       :: loop   !> TBD

    !! Local variables
    integer :: j, jj

    do j = 1,3
       if (kk(j) == +1) then
          khitp(j,loop) = khitp(j,loop) + 1
          do jj = 1,3
             vp(jj,j,loop) = vp(jj,j,loop) + dble(t(jj))
          end do
       else
          khite(j,loop) = khite(j,loop) + 1
          do jj = 1,3
             ve(jj,j,loop) = ve(jj,j,loop) + dble(t(jj))
          end do
       end if
    end do

  end subroutine accume

!! ================================================================

end module zeno_stat_utils

!!! ================================================================

!!! Local Variables:
!!! mode: f90
!!! time-stamp-line-limit: 30
!!! fill-column: 72
!!! End:
