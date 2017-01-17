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
!! This file provides the `zeno_problem_geometry` module which defines the
!! geometry of a program run as a whole.
!!
!! @author Mark Mansfield (formerly of Stevens Institute) and students
!! @author Walid Keyrouz (walid.keyrouz@nist.gov)
!! @date   Fri Nov  1 13:01:06 2013 EDT
!! 
!! @defgroup zeno_problem_geometry Problem Geometry for ZENO
!! 
!! Groups together routines that relate to the geometry of a problem to
!! be solved in a run of the program.
!! 
!! @todo Rewrite to use an object-oriented `shape` class hierarchy.
!! 
! Time-stamp: <2015-01-25 12:56:53 wtk>
! 
! ================================================================

module zeno_problem_geometry

  !! ================================================================

  use zeno_shape_distances, only : &
       & maxcube, procube, mincube, &
       & maxpillar, propillar, minpillar, &
       & maxsphere, prosphere, minsphere, &
       & maxtriangle, protriangle, mintriangle, &
       & maxdisk, prodisk, mindisk, &
       & maxcylinder, procylinder, mincylinder, &
       & maxtorus, protorus, mintorus, &
       & maxsolcyl, prosolcyl, minsolcyl, &
       & maxellipsoid, proellipsoid, minellipsoid

  implicit none
  private

  public do_launch

contains

  !! ================================================================
  !!
  !! Public routines
  !!
  !! ================================================================

  !> @brief Determine launch radius and enveloping box (original
  !! comment)
  !! 
  !! Compute the radius of the launch sphere and the enveloping box that
  !! are centered at the origin.
  !!
  !! @warning may need to be renamed to reflect functionality.
  !!
  !! @warning may need to move to a different module (e.g., problem
  !! geometry).
  !! 
  !! @warning too many parameters
  !! 
  !! @todo Choose a better algorithm and do not restrict the enclosing
  !! sphere to be centered at the origin.

  subroutine do_launch(maxelts, nelts, eltype, bv, rlaunch, &
       &       rotations, xyzlow, xyzhih)

    use zeno_codes_data

    !! Declare arguments
    integer, intent(in) :: maxelts         !> capacity of array
    integer, intent(in) :: nelts           !> # of elements processed
    integer, dimension(maxelts) :: eltype  !> type array of all elements
    real,    dimension(maxelts, 12) :: bv  !> TBD
    real                :: rlaunch         !> Radius of launch sphere
    real,    dimension(maxelts,3,3) :: rotations !> TBD
    real,    dimension(3) :: xyzlow              !> TBD
    real,    dimension(3) :: xyzhih              !> TBD

    !! Local variables
    real, dimension(3)   :: xaxis, yaxis, zaxis

    integer :: i
    real    :: dd
    real    :: sx1, sx2, sy1, sy2, sz1, sz2

    rlaunch = 0.0
    do i = 1,3
       xaxis(i) = 0.0
       yaxis(i) = 0.0
       zaxis(i) = 0.0
    end do
    xaxis(1) = 1.0
    yaxis(2) = 1.0
    zaxis(3) = 1.0

    Do_Nelts : do i = 1,nelts

       call do_launch_shape(eltype(i), i, maxelts, bv, rotations, &
            &               xaxis, yaxis, zaxis,                  &
            &               dd, sx1, sx2, sy1, sy2, sz1, sz2)

       rlaunch = amax1(dd,rlaunch)

       if (i == 1) then
          xyzlow(1) = sx1
          xyzhih(1) = sx2
          xyzlow(2) = sy1
          xyzhih(2) = sy2
          xyzlow(3) = sz1
          xyzhih(3) = sz2
       else
          xyzlow(1) = amin1(xyzlow(1),sx1)
          xyzhih(1) = amax1(xyzhih(1),sx2)
          xyzlow(2) = amin1(xyzlow(2),sy1)
          xyzhih(2) = amax1(xyzhih(2),sy2)
          xyzlow(3) = amin1(xyzlow(3),sz1)
          xyzhih(3) = amax1(xyzhih(3),sz2)
       end if

    end do Do_Nelts

  end subroutine do_launch

  !! ================================================================
  !!
  !! Private routine
  !!
  !! ================================================================

  !! Purpose of routine is to move horrendeous `if/switch` statement on
  !! `eltype(i)` out of loop over all or a subset of elements.
  !! 
  !! @warning Should rename to something more meaningful.
  !! 
  !! @warning too many parameters
  !! 
  !! @todo Use an object-oriented class hierachy

  subroutine do_launch_shape(elt_type, i, maxelts, bv, rotations, &
       &                     xaxis, yaxis, zaxis, &
       &                     dd, sx1, sx2, sy1, sy2, sz1, sz2)

    use zeno_codes_data

    !! Declare arguments
    integer, intent(in)  :: elt_type             !> Element type
    integer, intent(in)  :: i                    !> Element index
    integer, intent(in)  :: maxelts              !> capacity of array
    real,    dimension(maxelts, 12) :: bv        !> TBD
    real,    dimension(maxelts,3,3) :: rotations !> TBD
    real,    dimension(3):: xaxis, yaxis, zaxis
    real,    intent(out) :: dd, sx1, sx2, sy1, sy2, sz1, sz2

    !! Local variables
    real, dimension(3)   :: c, v1, v2, v3, n, n1, n2
    real, dimension(3)   :: x, y, zl, zh
    real, dimension(3,3) :: t
    integer :: j, k
    real    :: aa, al, bb, cc, r, r1, r2, side

    select case (elt_type)

    case (cube_code)
       do j = 1,3
          v1(j) = bv(i,j)
          side = bv(i,4)
       end do

       call maxcube(v1,side,dd)
       call procube(v1,side,xaxis,sx1,sx2)
       call procube(v1,side,yaxis,sy1,sy2)
       call procube(v1,side,zaxis,sz1,sz2)

    case (pillar_code)
       x(1) = bv(i,1)
       x(2) = bv(i,5)
       x(3) = bv(i,9)
       y(1) = bv(i,2)
       y(2) = bv(i,6)
       y(3) = bv(i,10)
       zl(1) = bv(i,3)
       zl(2) = bv(i,7)
       zl(3) = bv(i,11)
       zh(1) = bv(i,4)
       zh(2) = bv(i,8)
       zh(3) = bv(i,12)

       call maxpillar(x,y,zl,zh,dd)
       call propillar(x,y,zl,zh,xaxis,sx1,sx2)
       call propillar(x,y,zl,zh,yaxis,sy1,sy2)
       call propillar(x,y,zl,zh,zaxis,sz1,sz2)

    case (sphere_code)
       do j = 1,3
          c(j) = bv(i,j)
       end do
       r = bv(i,4)

       call maxsphere(c,r,dd)
       call prosphere(c,r,xaxis,sx1,sx2)
       call prosphere(c,r,yaxis,sy1,sy2)
       call prosphere(c,r,zaxis,sz1,sz2)

    case (triangle_code)
       do j = 1,3
          v1(j) = bv(i,j)
          v2(j) = bv(i,j+3)
          v3(j) = bv(i,j+6)
       end do

       call maxtriangle(v1,v2,v3,dd)
       call protriangle(v1,v2,v3,xaxis,sx1,sx2)
       call protriangle(v1,v2,v3,yaxis,sy1,sy2)
       call protriangle(v1,v2,v3,zaxis,sz1,sz2)

    case (disk_code)
       do j = 1,3
          c(j) = bv(i,j)
          n(j) = bv(i,j+3)
       end do
       r = bv(i,7)
       do j = 1,3
          do k = 1,3
             t(j,k) = rotations(i,j,k)
          end do
       end do

       call maxdisk(c,n,r,t,dd)
       call prodisk(c,n,r,t,xaxis,sx1,sx2)
       call prodisk(c,n,r,t,yaxis,sy1,sy2)
       call prodisk(c,n,r,t,zaxis,sz1,sz2)

    case (open_cylinder_code)
       do j = 1,3
          c(j) = bv(i,j)
          n(j) = bv(i,j+3)
       end do
       r = bv(i,7)
       al = bv(i,8)
       do j = 1,3
          do k = 1,3
             t(j,k) = rotations(i,j,k)
          end do
       end do

       call maxcylinder(c,n,r,al,t,dd)
       call procylinder(c,n,r,al,t,xaxis,sx1,sx2)
       call procylinder(c,n,r,al,t,yaxis,sy1,sy2)
       call procylinder(c,n,r,al,t,zaxis,sz1,sz2)

    case (solid_cylinder_code)
       do j = 1,3
          c(j) = bv(i,j)
          n(j) = bv(i,j+3)
       end do
       r = bv(i,7)
       al = bv(i,8)
       do j = 1,3
          do k = 1,3
             t(j,k) = rotations(i,j,k)
          end do
       end do

       call maxsolcyl(c,n,r,al,t,dd)
       call prosolcyl(c,n,r,al,t,xaxis,sx1,sx2)
       call prosolcyl(c,n,r,al,t,yaxis,sy1,sy2)
       call prosolcyl(c,n,r,al,t,zaxis,sz1,sz2)

    case (donut_code)
       do j = 1,3
          c(j) = bv(i,j)
          n(j) = bv(i,j+3)
       end do
       r1 = bv(i,7)
       r2 = bv(i,8)
       do j = 1,3
          do k = 1,3
             t(j,k) = rotations(i,j,k)
          end do
       end do

       call maxtorus(c,n,r1,r2,t,dd)
       call protorus(c,n,r1,r2,t,xaxis,sx1,sx2)
       call protorus(c,n,r1,r2,t,yaxis,sy1,sy2)
       call protorus(c,n,r1,r2,t,zaxis,sz1,sz2)

    case (ellipsoid_code)
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
             t(j,k) = rotations(i,j,k)
          end do
       end do

       call maxellipsoid(c,n1,n2,aa,bb,cc,t,dd)
       call proellipsoid(c,n1,n2,aa,bb,cc,t,xaxis,sx1,sx2)
       call proellipsoid(c,n1,n2,aa,bb,cc,t,yaxis,sy1,sy2)
       call proellipsoid(c,n1,n2,aa,bb,cc,t,zaxis,sz1,sz2)

    end select

  end subroutine do_launch_shape

  !! ================================================================

end module zeno_problem_geometry

!!! ================================================================

!!! Local Variables:
!!! mode: f90
!!! time-stamp-line-limit: 30
!!! fill-column: 72
!!! End:
