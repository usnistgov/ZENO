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
!! This file provides the `zeno_toss` module which defines routines that
!! perturb existing geometric entities.
!! 
!! @author Mark Mansfield (formerly of Stevens Institute) and students
!! @author Walid Keyrouz (walid.keyrouz@nist.gov)
!! @date   Mon Sep 30 15:58:33 2013 EDT
!! 
!! @defgroup zeno_toss Toss/Perturb Geometric Entities
!! 
!! Groups routines that perturb existing geometric entities.
!!
!! The module defines the following public entries: `toss_point`,
!! `toss_sphere`, `toss_cube`, `toss_pillar`, `toss_ellipsoid`,
!! `toss_cylinder`, and `toss_torus`
!!
!! @warning
!! The module's name may change to better reflect its functionality.
!! Also, the following exported routines are not implemented yet: 
!! `toss_pillar`, `toss_ellipsoid`, `toss_cylinder`, and `toss_torus`.
!!
!! @warning
!! Added hacks to get around compiler warnings about unsued dummy
!! arguments in the "not implemented yet" routines.
!!
!! @warning
!! Ideally, there should be a class hierarchy which defines a `toss`
!! method that then gets specialized by subclasses.
!! 
! Time-stamp: <2015-01-26 13:28:42 wtk>
!
! ================================================================

module zeno_toss

  !! ================================================================

  use zeno_rng, only : get_rand

  implicit none
  private

  public toss_point, toss_sphere, toss_cube, toss_pillar, &
       & toss_ellipsoid, toss_cylinder, toss_torus

contains

  !! ================================================================
  !!
  !! Public routines
  !!
  !! ================================================================

  !> To be filled in.

  subroutine toss_point(rt, maxelts, eltype, bv, nelts, &
       &                rotations, vaar, volume)

    use zeno_codes_data

    !! Declare arguments
    real, dimension(3)           :: rt         !> unknown yet!
    integer, intent(in)          :: maxelts    !> array capacity
    integer, dimension(maxelts)  :: eltype     !> elt. types array
    real, dimension(maxelts,12)  :: bv         !> unknown yet!
    integer, intent(in)          :: nelts      !> actual size
    real, dimension(maxelts,3,3) :: rotations  !> unknown yet!
    real, dimension(maxelts)     :: vaar       !> unknown yet!
    real                         :: volume     !> unknwon yet!

    !! Local variables
    real, dimension(3)   :: c, x, y, zl, zh, n, n1, n2
    real, dimension(3,3) :: tumble

    integer :: i, j, k, kdo
    real :: aa, bb, cc, alen, r, r1, r2, rad
    real :: side, sum, zip

    logical :: kdo_set

    zip = get_rand() * volume
    sum = 0.0

    kdo_set = .false.

    do i = 1,nelts
       sum = sum + vaar(i)
       if (zip < sum) then
          kdo = i
          kdo_set = .true.
          exit
       end if
    end do

    if (.not. kdo_set) then
       kdo = nelts
    end if

    if (eltype(kdo) == cube_code) then

       c(1) = bv(kdo,1)
       c(2) = bv(kdo,2)
       c(3) = bv(kdo,3)
       side = bv(kdo,4)
       call toss_cube(c,side,rt)

    else if (eltype(kdo) == pillar_code) then

       x(1) = bv(kdo,1)
       x(2) = bv(kdo,5)
       x(3) = bv(kdo,9)
       y(1) = bv(kdo,2)
       y(2) = bv(kdo,6)
       y(3) = bv(kdo,10)
       zl(1) = bv(kdo,3)
       zl(2) = bv(kdo,7)
       zl(3) = bv(kdo,11)
       zh(1) = bv(kdo,4)
       zh(2) = bv(kdo,8)
       zh(3) = bv(kdo,12)
       call toss_pillar(x,y,zl,zh,rt)

    else if (eltype(kdo) == sphere_code) then

       c(1) = bv(kdo,1)
       c(2) = bv(kdo,2)
       c(3) = bv(kdo,3)
       r = bv(kdo,4)
       call toss_sphere(c,r,rt)

    else if (eltype(kdo) == donut_code) then
        
       do j = 1,3
          c(j) = bv(i,j)
          n(j) = bv(i,j+3)
       end do
       r1 = bv(i,7)
       r2 = bv(i,8)
       do j = 1,3
          do k = 1,3
             tumble(j,k) = rotations(i,j,k)
          end do
       end do
       call toss_torus(c, n, r1, r2, tumble, rt)

    else if (eltype(kdo) == solid_cylinder_code) then

       do j = 1,3
          c(j) = bv(i,j)
          n(j) = bv(i,j+3)
       end do
       rad = bv(i,7)
       alen = bv(i,8)
       do j = 1,3
          do k = 1,3
             tumble(j,k) = rotations(i,j,k)
          end do
       end do
       call toss_cylinder(c,n,rad,alen,tumble,rt)

    else if (eltype(kdo) == ellipsoid_code) then

       do j = 1,3
          c(j) = bv(i,j)
          n1(j) = bv(i,j+3)
          n2(j) = bv(i,j+6)
       end do
       aa = bv(i,10)
       bb = bv(i,11)
       cc = bv(i,12)
       do j = 1,3
          do k = 1,3
             tumble(j,k) = rotations(i,j,k)
          end do
       end do
       call toss_ellipsoid(c,n1,n2,aa,bb,cc,tumble,rt)

    else
        
       stop 'bad code in toss_point'

    end if

  end subroutine toss_point

  !! ----------------------------------------------------------------

  !> To be filled in!

  subroutine toss_sphere(c, r, rt)

    !! Declare arguments
    real, dimension(3) :: c, rt
    real               :: r

    !! Local variables
    integer :: i
    real    :: sum

    do
       do i = 1,3
          rt(i) = (2.0*get_rand() - 1.0)*r
       end do
       sum = rt(1)**2 + rt(2)**2 + rt(3)**2
       if (sum <= r**2) exit
    end do

    do i = 1,3
       rt(I) = rt(i) + c(i)
    end do

  end subroutine toss_sphere

  !! ----------------------------------------------------------------

  !> To be filled in.

  subroutine toss_cube(c, side, rt)

    !! Declare arguments
    real, dimension(3) :: c, rt
    real               :: side

    !! Local variables
    integer :: i

    do i = 1,3
       rt(i) = get_rand()*side + c(i)
    end do

  end subroutine toss_cube

  !! ----------------------------------------------------------------

  !> To be filled in.
  !!  
  !! @warning Incomplete implementation!

  subroutine toss_pillar(x, y, zl, zh, rt)

    !!GCC$ ATTRIBUTES unused :: x, y, zl, zh, rt

    !! Declare arguments
    real, dimension(3) :: x, y, zl, zh, rt

    !! Get around compiler warning for unused arguments
    real, dimension(3) :: x_2, y_2, zl_2, zh_2, rt_2

    !! HACK !!
    if (.false.) then
       x_2(1) = x(1)
       y_2(2) = y(2)
       zl_2(3) = zl(3)
       zh_2(1) = zh(1)
       rt_2(2) = rt(2)
    end if

    write(*,"('SORRY!!!!   toss_pillar not ready')")

    stop

  end subroutine toss_pillar

  !! ----------------------------------------------------------------

  !> To be filled in.
  !! 
  !! @warning Incomplete implementation!

  subroutine toss_ellipsoid(c, n1, n2, aa, bb, cc, tumble, rt)

    !!GCC$ ATTRIBUTES unused :: c, n1, n2, aa, bb, cc, tumble, rt

    !! Declare arguments
    real                 :: aa, bb, cc
    real, dimension(3)   :: c, n1, n2, rt
    real, dimension(3,3) :: tumble

    !! Get around compiler warning for unused arguments
    real                 :: aa_2, bb_2, cc_2
    real, dimension(3)   :: c_2, n1_2, n2_2, rt_2
    real, dimension(3,3) :: tumble_2

    !! HACK !!
    if (.false.) then
       aa_2 = aa
       bb_2 = bb
       cc_2 = cc
       c_2(1) = c(1)
       n1_2(2) = n1(2)
       n2_2(3) = n2(3)
       rt_2(1) = rt(1)
       tumble_2(1,1) = tumble(1,1)
    end if

    write(*,"('SORRY  ---  toss_ellipsoid not ready')")

    stop
  end subroutine toss_ellipsoid

  !! ----------------------------------------------------------------

  !> To be filled in.
  !! 
  !! @warning Incomplete implementation!

  subroutine toss_cylinder(c, n, rad, alen, tumble, rt)

    !!GCC$ ATTRIBUTES unused :: c, n, rad, alen, tumble, rt)

    !! Declare arguments
    real                 :: rad, alen
    real, dimension(3)   :: c, n, rt
    real, dimension(3,3) :: tumble

    !! Declare arguments
    real                 :: rad_2, alen_2
    real, dimension(3)   :: c_2, n_2, rt_2
    real, dimension(3,3) :: tumble_2

    !! HACK !!
    if (.false.) then
       rad_2 = rad
       alen_2 = alen
       c_2(1) = c(1)
       n_2(2) = n(2)
       rt_2(3) = rt(3)
       tumble_2(1,1) = tumble(1,1)
    end if

    write(*,"('SORRY  ---  toss_cylinder not ready')")

    stop
  end subroutine toss_cylinder

  !! ----------------------------------------------------------------

  !> To be filled in.
  !! 
  !! @warning Incomplete implementation!

  subroutine toss_torus(c, n, r1, r2, tumble, rt)

    !!GCC$ ATTRIBUTES unused :: n, r1, r2, tumble, rt

    !! Declare arguments
    real                 :: r1, r2
    real, dimension(3)   :: c, n, rt
    real, dimension(3,3) :: tumble

    !! Get around compiler warning for unused arguments
    real                 :: r1_2, r2_2
    real, dimension(3)   :: c_2, n_2, rt_2
    real, dimension(3,3) :: tumble_2

    !! HACK !!
    if (.false.) then
       r1_2 = r1
       r2_2 = r2
       c_2(1) = c(1)
       n_2(2) = n(2)
       rt_2(3) = rt(3)
       tumble_2(1,1) = tumble(1,1)
    end if

    write(*,"('SORRY!!!!   toss_torus not ready')")

    stop
  end subroutine toss_torus

  !! ================================================================

end module zeno_toss

!!! ================================================================

!!! Local Variables:
!!! mode: f90
!!! time-stamp-line-limit: 30
!!! fill-column: 72
!!! End:
