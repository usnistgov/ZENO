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
!! This file defines the `zeno_sort_utils` module which groups some
!! sorting routines.
!!
!! @author Mark Mansfield (formerly of Stevens Institute) and students
!! @author Walid Keyrouz <walid.keyrouz@nist.gov>
!! @date   Fri Oct 25 10:34:03 2013 EDT
!! 
!! @defgroup zeno_sort_utils Sorting Utilities
!!
!! Groups together sort utility routines.
!!
! Time-stamp: <2015-02-23 18:07:05 wtk>
! 
! ================================================================

module zeno_sort_utils

  !! ================================================================

  implicit none

  private
  public sort, sort3, listersort

contains

  ! ================================================================
  !
  ! Public routines
  !
  ! ================================================================

  !> @brief Computes the squares of a 3-element array and sorts it
  !! using `bubble-sort`!

  subroutine sort(a, a2s, ndif)

    !! Declare arguments
    real, dimension(3)   :: a, a2s
    integer, intent(out) :: ndif

    !! Local variables
    integer :: i, nwk
    real :: z

    do i = 1,3
       a2s(i) = a(i)**2
    end do

    ndif = 3
    nwk = 1

    do while (nwk < ndif)
       if (nwk == 0) nwk = 1
       if (a2s(nwk) < a2s(nwk+1)) then
          nwk = nwk + 1
       else if (a2s(nwk) > a2s(nwk+1)) then
          z = a2s(nwk)
          a2s(nwk) = a2s(nwk+1)
          a2s(nwk+1) = z      
          nwk = nwk - 1
       else
          a2s(nwk+1) = a2s(ndif)
          ndif = ndif - 1
       end if
    end do

  end subroutine sort

  !! ----------------------------------------------------------------

  !> Sort three parameters s.t. `as > bs > cs` using `bubble-sort`.

  subroutine sort3(as, bs, cs)

    !! Declare arguments
    real, intent(in out) :: as, bs, cs    !> Scalars to sort

    !! Local variables
    real, dimension(3) :: a
    integer :: nwk
    real    :: save

    a(1) = as
    a(2) = bs
    a(3) = cs

    nwk = 1
    do while (nwk < 3)
       if (nwk == 0) nwk = 1
       if (a(nwk) >= a(nwk+1)) then
          nwk = nwk + 1
       else
          save = a(nwk)
          a(nwk) = a(nwk+1)
          a(nwk+1) = save
          nwk = nwk - 1
       end if
    end do

    as = a(1)
    bs = a(2)
    cs = a(3)

  end subroutine sort3

  !! ----------------------------------------------------------------

  !> @brief Seems to sort `rlist` and `nebtab` via `bubblesort`!
  !!
  !! @warning should have a more meaningful name!
  !!
  !! @warning Routine uses @c bubblesort for sorting.  This may be a
  !! performance bottleneck.  Must check the values of @c ninn and @c
  !! nneb which are the sizes of @c rlist and @c nebtab.
  !!
  !! Mitigating factor re using @c bubblesort is that the number of
  !! elements being sorted, @c ninn, is always <= 4.

  subroutine listersort(rlist, ninn, nneb, nebtab, maxelts)

    use zeno_debug

    !! Declare arguments
    real, dimension(maxelts)    :: rlist   !> TBD
    integer, intent(in out)     :: ninn    !> Actual num. of neighbors
    integer, intent(in)         :: nneb    !> Max num. of neighbors
    integer, dimension(maxelts) :: nebtab  !> Neighbors list
    integer, intent(in)         :: maxelts !> Num. of elements

    !! Local variables
    integer :: kopy, nwk
    real    :: copy

    !> Update @c log data
    calls_LS_cnt = calls_LS_cnt + 1

    nwk = 1
    do while (nwk < ninn)
       if (nwk == 0) nwk = 1
       if (rlist(nwk) <= rlist(nwk+1)) then
          nwk = nwk + 1
       else
          copy = rlist(nwk)
          rlist(nwk) = rlist(nwk+1)
          rlist(nwk+1) = copy
          kopy = nebtab(nwk)
          nebtab(nwk) = nebtab(nwk+1)
          nebtab(nwk+1) = kopy
          nwk = nwk - 1
       end if
    end do

    if (ninn > nneb) ninn = nneb

  end subroutine listersort

  !! ================================================================

end module zeno_sort_utils

!!! Local Variables:
!!! time-stamp-line-limit: 30
!!! mode: f90
!!! End:
