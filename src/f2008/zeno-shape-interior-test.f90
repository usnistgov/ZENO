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
!! This file provides the `zeno_shape_interior_test` module which defines
!! routines for determining if a point is inside a particular shape.
!! 
!! @warning Need definitions for all shapes
!! 
!! @todo Code crying out for an object-oriented aproach with a `shape`
!! class hierarchy.
!!
!! @author Mark Mansfield (formerly of Stevens Institute) and students
!! @author Walid Keyrouz (walid.keyrouz@nist.gov)
!! @date   Thu Oct 24 16:40:48 2013 EDT
!! 
!! @defgroup zeno_shape_interior_test Shape Interior Test
!! 
!! Groups routines that determine if a point is inside a particular
!! shape.
!!
!! @warning Need definitions for all shapes.
!! 
!! @warning Should be merged as a method in a `shape` class hierarchy.
!!
!! @todo Use an object-oriented approach with a class hierarchy.
!!
! Time-stamp: <2015-01-05 16:45:00 walid>
! 
! ================================================================

module zeno_shape_interior_test

  !! ================================================================

  use zeno_vectors, only: vector_difference, rotate, pythag

  implicit none
  private

  public inbody, inside_shape, &
       & insidecube, insidepillar, insidesphere, insidetorus, &
       & insidesolcyl, insideellipsoid

contains

  !! ================================================================
  !!
  !! Public routines
  !!
  !! ================================================================

  !> Returns inside as true if the point `rt` lies inside the body
  !!
  !! @warning Routine does not handle all shape types!

  subroutine inbody(rt, maxelts, eltype, bv, nelts, rotations, &
       &            inside, early, mess)

    use zeno_codes_data

    !! Declare arguments
    real, dimension(3)           :: rt        !> Point being tested
    integer, intent(in)          :: maxelts   !> Array capacity
    integer, dimension(maxelts)  :: eltype    !> Elements array
    real, dimension(maxelts,12)  :: bv        !> TBD
    integer, intent(in)          :: nelts     !> Number of elements
    real, dimension(maxelts,3,3) :: rotations !> Element rotations
    logical, intent(out)         :: inside    !> test result
    logical, intent(out)         :: early     !> TBD
    character(len=10)            :: mess      !> Message?

    !! Local variables
    integer :: i

    inside = .false.
    early = .false.

    !> WK: should not have to iterate over all elements with a @c
    !>     proper spatial data structure.

    do_nelts: do i = 1,nelts

       call inside_shape(rt, eltype(i), i, maxelts, bv, rotations, &
            &            inside, early, mess)
       if (inside .or. early) return

    end do do_nelts

  end subroutine inbody

  !! ----------------------------------------------------------------

  !> Tests if a point `p` is inside a particular shape.  Used to get
  !! rid of horrendous `if` statement on the type of an element.

  subroutine inside_shape(rt, elt_type, i, maxelts, bv, rotations, &
       &                  inside, early, mess)

    use zeno_codes_data

    !! Declare arguments
    real, dimension(3)  :: rt        !> Point being tested
    integer, intent(in) :: elt_type
    integer, intent(in) :: i
    integer, intent(in) :: maxelts
    real, dimension(maxelts,12)  :: bv
    real, dimension(maxelts,3,3) :: rotations !> Element rotations
    logical, intent(out) :: inside
    logical, intent(out) :: early
    character(len=10)    :: mess

    !! Local variables
    real, dimension(3)   :: x, y, zl, zh
    real, dimension(3)   :: c, n, n1, n2
    real, dimension(3,3) :: t
    integer :: j, k
    real    :: aa, bb, cc, alen
    real    :: r, r1, r2, rad, side

    inside = .false.
    early  = .false.

    select case (elt_type)

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
       call insidepillar(x,y,zl,zh,rt,inside)

    case (cube_code)
       do j = 1,3
          c(j) = bv(i,j)
       end do
       side = bv(i,4)
       call insidecube(c,side,rt,inside)

    case (sphere_code)
       do j = 1,3
          c(j) = bv(i,j)
       end do
       r = bv(i,4)
       call insidesphere(c,r,rt,inside)

    case (disk_code)
       mess = 'DISK      '
       early = .true.

    case (open_cylinder_code)
       mess = 'O-CYLINDER'
       early = .true.

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
       call insidetorus(c,n,r1,r2,t,rt,inside)

    case (solid_cylinder_code)
       do j = 1,3
          c(j) = bv(i,j)
          n(j) = bv(i,j+3)
       end do
       rad = bv(i,7)
       alen = bv(i,8)
       do j = 1,3
          do k = 1,3
             t(j,k) = rotations(i,j,k)
          end do
       end do
       call insidesolcyl(c,n,rad,alen,t,rt,inside)

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
       call insideellipsoid(c,n1,n2,aa,bb,cc,t,rt,inside)

    end select

  end subroutine inside_shape

  !! ----------------------------------------------------------------

  !> @brief Tests if point `p` is inside a cube.
  !!
  !! Given a cube with low corner at `v1` and with length of side =
  !! `side`, also given an arbitrary point `p` report in 
  !! `"inside"` whether or not `p` is inside the cube

  subroutine insidecube(v1, side, p, inside)

    !! Declare arguments
    real, dimension(3)   :: v1        !> Lower corner of cube
    real, intent(in)     :: side      !> Size of cube
    real, dimension(3)   :: p         !> Point tested
    logical, intent(out) :: inside    !> Test result

    !! Local variables
    integer :: i, nwk
    real    :: s2, ss
    real, dimension(3) :: c, pr

    s2 = side/2.0

    do i = 1,3
       c(i) = v1(i) + s2
    end do

    do i = 1,3
       pr(i) = p(i) - c(i)
    end do

    do i = 1,3
       pr(i) = abs(pr(i))
    end do

    nwk = 1
    do while (nwk < 3)
       if (pr(nwk) <= pr(nwk+1)) then
          nwk = nwk + 1
       else
          ss = pr(nwk)
          pr(nwk) = pr(nwk+1)
          pr(nwk+1) = ss
          nwk = nwk - 1
       end if
       if (nwk == 0) nwk = 1
    end do

    inside = pr(3) < s2

  end subroutine insidecube

  !! ----------------------------------------------------------------

  !> @brief Tests if point `p` is inside the `pillar` body.
  !!
  !! Is the point [p(1),p(2)] inside the triangle [x(1),y(1)],
  !! [x(2),y(2)], [x(3),y(3)]?
  !!
  !! @warning Description from original comment in code seems to be
  !! misleading!

  subroutine insidepillar(x, y, zl, zh, p, inside)

    !! Declare arguments
    real, dimension(3) :: x, y, zl, zh    !> Pillar specs
    real, dimension(3) :: p               !> Point tested
    logical, intent(out) :: inside        !> Test result

    !! Local variables
    real, dimension(2) :: b1, b2, s, alpha
    real               :: zlower, zoom, zupper

    b1(1) = x(2) - x(1)
    b1(2) = y(2) - y(1)
    b2(1) = x(3) - x(1)
    b2(2) = y(3) - y(1)
    s(1) = p(1) - x(1)
    s(2) = p(2) - y(1)
    zoom = b1(1)*b2(2) - b2(1)*b1(2)
    alpha(1) = (s(1)*b2(2) - b2(1)*s(2))/zoom
    alpha(2) = (b1(1)*s(2) - s(1)*b1(2))/zoom

    if (alpha(1) < 0.0) then
       inside = .false.
       return
    end if

    if (alpha(2) < 0.0) then
       inside = .false.
       return
    end if

    if (alpha(1)+alpha(2) > 1.0) then
       inside = .false.
       return
    end if

    !> If we have fallen out here, then we know that the image of the
    !> point falls inside the triangle.
        
    zupper = zh(1) + alpha(1)*(zh(2)-zh(1)) + alpha(2)*(zh(3)-zh(1))
    zlower = zl(1) + alpha(1)*(zl(2)-zl(1)) + alpha(2)*(zl(3)-zl(1))

    inside = (zlower <= p(3)) .and. (p(3) <= zupper)

  end subroutine insidepillar

  !! ----------------------------------------------------------------

  !> @brief tests if a point `p` is inside a sphere.
  !!
  !! Given a sphere (center `c`, radius `r`) and an arbitrary point
  !! `p`, return `inside=.true.` if `p` is inside the sphere.

  subroutine insidesphere(c, r, p, inside)

    !! Declare arguments
    real, dimension(3)   :: c         !> Center of sphere
    real, intent(in)     :: r         !> Radius of sphere
    real, dimension(3)   :: p         !> Point tested
    logical, intent(out) :: inside    !> Test result

    !! Local variables
    real :: r1

    call pythag(c,p,r1)
    inside = r1 < r

  end subroutine insidesphere

  !! ----------------------------------------------------------------

  !> @brief Tests if point `p` is inside a torus.
  !!
  !! @verbatim
  !! TORUS[(cx,cy,cz),(nx,ny,nz),r1,r2]
  !! 
  !! Given a torus (center `c`, unit normal `n`, radii `r1` and 
  !! `r2`), return `inside` as `true` if `p` is inside the torus.
  !! 
  !!           x                         x
  !!       x       x                 x       x
  !!      x         x               x         x   ___
  !!      x         x               x         x    |
  !!       x       x                 x       x     |  r2
  !!           x                         x        ---
  !! 
  !!           |--------2 * r1 ----------|
  !! @endverbatim

  subroutine insidetorus(c, n, r1, r2, t, p, inside)

    !!GCC$ ATTRIBUTES unused :: n

    !! Declare arguments
    real, dimension(3)   :: c, n     !> Torus center and unit normal
    real, intent(in)     :: r1, r2   !> Torus radii
    real, dimension(3)   :: p        !> Point tested
    logical, intent(out) :: inside   !> Test result

    !! Local variables
    real, dimension(3)   :: pmc, q
    real, dimension(3,3) :: t
    real                 :: d1, qq

    !! Get around compiler warning for unused dummy argument(s)
    real, dimension(3) :: n_2

    !! HACK !!
    if (.false.) then
       n_2(1) = n(1)
    end if

    call vector_difference(p,c,pmc)
    call rotate(q,t,pmc)

    qq = sqrt(q(1)**2+q(2)**2)

    d1 = (qq-r1)**2 + q(3)**2
    d1 = sqrt(d1)
    if (d1 < r2) then
       inside = .true.
       return
    end if

    d1 = (qq+r1)**2 + q(3)**2
    d1 = sqrt(d1)
    if (d1 < r2) then
       inside = .true.
       return
    end if

    inside = .false.

  end subroutine insidetorus

  !! ----------------------------------------------------------------

  !> @brief tests if a point `p` is inside a solid cylinder.
  !!
  !! @warning Misleading description comment in original code.
  !!
  !! Given a solid cylinder (center `c`, unit normal `n`, radius 
  !! `r`, length `l`), returns the minimum distance between points on
  !! the cylinder and the arbitrary point `p`.
  !!
  !! @warning There is no `insidecylinder` routine!
 
  subroutine insidesolcyl(c, n, r, l, t, p, inside)

    !!GCC$ ATTRIBUTES unused :: n

    !! Declare arguments
    real, dimension(3)   :: c, n      !> Cylinder center and unit normal
    real, intent(in)     :: r, l      !> Cylinder radius and length
    real, dimension(3)   :: p         !> Point being tested
    logical, intent(out) :: inside    !> Test result

    !! Local variables
    real, dimension(3,3) :: t
    real, dimension(3)   :: pmc, q
    real                 :: dp, ro

    !! Get around compiler warning for unused dummy argument(s)
    real, dimension(3) :: n_2

    !! HACK !!
    if (.false.) then
       n_2(1) = n(1)
    end if

    call vector_difference(p,c,pmc)
    call rotate(q,t,pmc)
    ro = sqrt(q(1)**2 + q(2)**2)
    dp = abs(ro-r)

    inside = .false.

    if (abs(q(3)) < l/2.0) then
       if (ro <= r) then
          inside = .true.
       end if
    end if

  end subroutine insidesolcyl

  !! ----------------------------------------------------------------

  !> @brief tests if point `p` is inside ellipsoid.
  !!
  !! Return `inside=.true.` if point `p` is on the interior of the
  !! ellipsoid

  subroutine insideellipsoid(c, n1, n2, aa, bb, cc, t, po, inside)

    !!GCC$ ATTRIBUTES unused :: n1, n2

    !! Declare arguments
    real, dimension(3)   :: c, n1, n2  !> Ellipsoid center and two normals
    real, intent(in)     :: aa, bb, cc !> Ellipsoid axial lengths
    real, dimension(3,3) :: t          !> Special rotation matrix
    real, dimension(3)   :: po         !> Point being tested
    logical, intent(out) :: inside     !> Test result

    !! Local variables
    real, dimension(3) :: pomc, p
    real               :: zz

    !! Get around compiler warning for unused dummy argument(s)
    real, dimension(3) :: n1_2, n2_2

    !! HACK !!
    if (.false.) then
       n1_2(1) = n1(1)
       n2_2(2) = n2(2)
    end if

    call vector_difference(po,c,pomc)
    call rotate(p,t,pomc)

    zz = (p(1)/aa)**2 + (p(2)/bb)**2 + (p(3)/cc)**2

    inside = zz < 1.0

  end subroutine insideellipsoid

  !! ================================================================

end module zeno_shape_interior_test

!!! ================================================================

!!! Local Variables:
!!! mode: f90
!!! time-stamp-line-limit: 30
!!! fill-column: 72
!!! End:
