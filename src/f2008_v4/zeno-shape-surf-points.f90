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
!! This file provides the `zeno_shape_surface_points` module which defines
!! routines that generate a uniform point coverage for shape surfaces.
!! 
!! @warning Need definitions for all shapes
!! 
!! @todo Code crying out for an object-oriented aproach with a `shape`
!! class hierarchy.
!!
!! @author Mark Mansfield (formerly of Stevens Institute) and students
!! @author Walid Keyrouz (walid.keyrouz@nist.gov)
!! @date   Fri Oct 25 11:39:21 2013 EDT
!! 
!! @defgroup zeno_shape_surface_points Shape Surface Points
!! 
!! Groups routines that generate uniformly distributed points on shape
!! surfaces.
!! 
!! @warning should be merged as methods in a `shape` class hierarchy.
!! 
! Time-stamp: <2015-01-26 12:30:41 wtk>
! 
! ================================================================

module zeno_shape_surface_point

  !! ================================================================

  use numeric_kinds
  use numeric_constants

  use zeno_sphere, only : sphere
  use zeno_vectors, only : &
       & vector_difference, dotproduct, pythag0, makepolar, backtransform
  use zeno_rng, only : get_rand, jrand
  use zeno_shape_interior_test, only: inside_shape
  use zeno_shape_surface_area, only : pillarsurf
  use zeno_sort_utils, only : sort3

  implicit none
  private

  public getsurface

contains

  !! ================================================================
  !!
  !! Public routines
  !!
  !! ================================================================

  !> Generate a point, `p1`, distributed uniformly over the surface
  !!
  !! @warning may not handle all shapes.

  subroutine getsurface(maxelts, eltype, bv, nelts, saar, total, p1, &
    &                   trials, rotations, loop)

    use zeno_codes_data

    !! Declare arguments
    integer, intent(in)          :: maxelts   !> Array capacity
    integer, dimension(maxelts)  :: eltype    !> Elements array
    real, dimension(maxelts,12)  :: bv        !> TBD
    integer, intent(in)          :: nelts     !> Array size
    real, dimension(maxelts)     :: saar      !> TBD
    real, intent(out)            :: total     !> TBD
    real, dimension(3)           :: p1        !> Generated point
    real(dp_k), dimension(20)    :: trials    !> TBD
    real, dimension(maxelts,3,3) :: rotations !> TBD
    integer, intent(in out)      :: loop      !> TBD

    !! Local variables
    logical :: inside = .false.
    integer :: i, kdo
    real    :: sum, zip

    logical :: kdo_set

    logical :: early
    character(len=10) :: mess

    Rpt_While_Outside: do

       trials(loop) = trials(loop) + 1.0d0

       zip = get_rand()*total
       sum = 0.0

       kdo_set = .false.
       do i = 1,nelts
          sum = sum + saar(i)
          if (zip < sum) then
             kdo = i; kdo_set = .true.
             exit
          end if
       end do

       if (.not. kdo_set) kdo = nelts

       i = kdo

       call dover_shape(eltype(kdo), i, maxelts, bv, p1)

       !> p1 is now a point distributed over the surface
       !> of the elements.  But we must remove it if it
       !> is inside any other element.

       !> WK: should not have to iterate over all elements with a
       !>     proper spatial data structure

       Do_Nelts_I: do i = 1,nelts

          if (i /= kdo) then

             call inside_shape(p1, eltype(i), i, maxelts, bv, rotations, &
                  &            inside, early, mess)

             if (inside) exit

          end if                   ! if (i /= kdo)

       end do Do_Nelts_I

       if (.not. inside) exit

    end do Rpt_While_Outside

  end subroutine getsurface

  !! ================================================================
  !! 
  !! Private routines
  !! 
  !! ================================================================

  !> Generate point distributed over a particular shape surf.
  !! Purpose of routine is to move @e horrendous `if/switch`
  !! statement on `eltype` out of loop over all or a subset of
  !! elements.
  !! 
  !! @warning should be renamed to something more meaninful.

  subroutine dover_shape(elt_type, i, maxelts, bv, p)

    use zeno_codes_data

    !! Declare arguments
    integer, intent(in)   :: elt_type
    integer, intent(in)   :: i
    integer, intent(in)   :: maxelts
    real, dimension(maxelts,12) :: bv
    real, dimension(maxelts,3,3) :: rotations !> TBD
    real, dimension(3)    :: p

    !! Local variables
    real, dimension(3) :: c, n, v1, v2, v3
    real, dimension(3) :: x, y, zl, zh
    real, dimension(3) :: n1, n2
    real, dimension(3,3) :: t

    integer :: j, k
    real    :: aa, al, bb, cc, r, r1, r2
    real    :: side

    select case (elt_type)
       
    case (cube_code)
       do j = 1,3
          c(j) = bv(i,j)
       end do
       side = bv(i,4)
       call dovercube(c,side,p)

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
       call doverpillar(x,y,zl,zh,p)

    case (sphere_code)
       do j = 1,3
          c(j) = bv(i,j)
       end do
       r = bv(i,4)
       call doversphere(c,r,p)

    case (triangle_code)
       do j = 1,3
          v1(j) = bv(i,j)
          v2(j) = bv(i,j+3)
          v3(j) = bv(i,j+6)
       end do
       call dovertriangle(v1,v2,v3,p)

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
       call doverdisk(c,n,r,t,p)

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
       call dovercylinder(c,n,r,al,t,p)

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
       call doversolcyl(c,n,r,al,t,p)

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
       call dovertorus(c,n,r1,r2,t,p)

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
       call doverellipsoid(c,n1,n2,aa,bb,cc,t,p)

    end select

  end subroutine dover_shape

  !! ----------------------------------------------------------------

  !> Generate point distributed over a cube.
  !! 
  !! @warning should be renamed to something more meaningful.

  subroutine dovercube(c, side, p1)

    !! Declare arguments
    real, dimension(3) :: c     !> Cube lower corner
    real, intent(in)   :: side  !> Cube side length
    real, dimension(3) :: p1    !> Point generated

    !! Local variables
    integer :: idir, is, kk

    idir = jrand(1,3)

    do kk = 1,3
       if (kk /= idir) then
          p1(kk) = side*(get_rand()-0.5)
       else
          is = jrand(1,2)
          if (is == 1) then
             p1(kk) = -side/2.0
          else
             p1(kk) = side/2.0
          end if
       end if
       p1(kk) = p1(kk) + c(kk) + side/2.0
    end do

  end subroutine dovercube

  !! ----------------------------------------------------------------

  !> Generate a point distributed randomly over sphere
  !! 
  !! @warning should be renamed!

  subroutine doversphere(c, r, p1)

    !! Declare arguments
    real, dimension(3) :: c     !> Sphere center
    real, intent(in)   :: r     !> Sphere radius
    real, dimension(3) :: p1    !> Point generated

    !! Local variables
    integer :: i

    call sphere(p1,r)
    do i = 1,3
       p1(i) = p1(i) + c(i)
    end do

  end subroutine doversphere

  !! ----------------------------------------------------------------

  !> Generate a point distributed randomly over a triangle
  !>
  !> @warning should be renamed.
  !>
  !> @warning spaghetti code

  subroutine dovertriangle(p1, p2, p3, tt)

    !! Declare arguments
    real, dimension(3) :: p1, p2, p3    !> Triangle vertices
    real, dimension(3) :: tt            !> Generated point

    !! Local variables
    real, dimension(3) :: p23, p21, q, av, p31
    real, dimension(3) :: t1
    integer :: i
    real    :: a, aa1, aa2, bb1, bb2, bot, cc1, cc2
    real    :: gamma2, gamma3, randy
    real    :: tat, th, tl, top, ts

    call vector_difference(p2,p3,p23)   ! nm
    call vector_difference(p2,p1,p21)   ! nm
    call dotproduct(p23,p21,top)        ! nm**2
    call dotproduct(p21,p21,bot)        ! nm**2
    tat = top/bot                       ! nm**0
    do i = 1,3
       q(i) = tat*p1(i) + (1.0-tat)*p2(i)  ! nm
    end do

    !> av = altitude vector normal to (p1...p2) side
    call vector_difference(p3,q,av)   ! nm
    !> a = altitude normal to (p1...p2) side
    call pythag0(av,a)             ! nm

    th = amax1(1.0,tat)         ! nm**0
    tl = amin1(0.0,tat)         ! nm**0

    do
       ts = tl + get_rand()*(th-tl)        ! nm**0
       randy = get_rand()
       do i = 1,3
          tt(i) = ts*p1(i) + (1.0-ts)*p2(i) + randy*av(i)   ! nm
       end do

       call vector_difference(tt,p1,t1)        ! nm
       call dotproduct(t1,p21,aa1)             ! nm**2
       bb1 = bot                               ! nm**2
       call vector_difference(p3,p1,p31)       ! nm
       call dotproduct(p31,p21,cc1)            ! nm**2
       call dotproduct(t1,p31,aa2)             ! nm**2
       bb2 = cc1                               ! nm**2
       call dotproduct(p31,p31,cc2)            ! nm**2

       gamma2 = (aa1*cc2 - cc1*aa2)/(bb1*cc2-cc1*bb2)   ! nm**0
       if (gamma2 < 0.0) cycle

       gamma3 = (bb1*aa2 - aa1*bb2)/(bb1*cc2-cc1*bb2)   ! nm**0
       if (gamma3 < 0.0) cycle
       if (gamma2+gamma3 > 1.0) cycle

       exit
    end do

  end subroutine dovertriangle

  !! ----------------------------------------------------------------

  !> Generate point distributed randomly over circular disk
  !! 
  !! @warning should be renamed
  !! 
  !! @warning spaghetti code

  subroutine doverdisk(c, n, r, t, p1)

    !!GCC$ ATTRIBUTES unused :: n

    !! Declare arguments
    real, dimension(3)   :: c, n    !> Disk center and unit normal
    real, intent(in)     :: r       !> Disk radius
    real, dimension(3,3) :: t       !> Special rotation matrix
    real, dimension(3)   :: p1      !> Generated point

    !! Local variables
    real :: xx

    !! Get around compiler warning for unused dummy argument(s)
    real, dimension(3) :: n_2

    !! HACK !!
    if (.false.) then
       n_2(1) = n(1)
    end if

    do
       p1(1) = 2.0*r*(get_rand()-0.5)
       p1(2) = 2.0*r*(get_rand()-0.5)
       xx = p1(1)**2 + p1(2)**2
       xx = sqrt(xx)
       if (xx <= r) exit
    end do

    p1(3) = 0.0

    call backtransform(c,t,p1)

  end subroutine doverdisk

  !! ----------------------------------------------------------------

  !> Distribute a point randomly over cylinder
  !! 
  !! @warning should be renamed
  !! 
  !! @warning uses anonymous numeric constants

  subroutine dovercylinder(c, n, r, al, t, p1)

    !!GCC$ ATTRIBUTES unused :: n

    !! Declare arguments
    real, dimension(3)   :: c, n    !> Cylinder center and unit normal
    real, intent(in)     :: r, al   !> Cylinder radius and length
    real, dimension(3,3) :: t       !> Special rotation matrix
    real, dimension(3)   :: p1      !> Generated point

    !! Local variables
    real :: theta

    !! Get around compiler warning for unused dummy argument(s)
    real, dimension(3) :: n_2

    !! HACK !!
    if (.false.) then
       n_2(1) = n(1)
    end if

    p1(3) = al * (get_rand()-0.5)
    theta = 2.0 * M_PI_sp * get_rand()
    p1(1) = r * cos(theta)
    p1(2) = r * sin(theta)

    call backtransform(c,t,p1)

  end subroutine dovercylinder

  !! ----------------------------------------------------------------

  !> Distribute a point randomly over a solid cylinder
  !! 
  !! @warning should be renamed
  !! 
  !! @warning uses anonymous numeric constants

  subroutine doversolcyl(c, n, r, al, t, p1)

    !! Declare arguments
    real, dimension(3)   :: c, n    !> Cylinder center and unit normal
    real, intent(in)     :: r, al   !> Cylinder radius and length
    real, dimension(3,3) :: t       !> Special rotation matrix
    real, dimension(3)   :: p1      !> Generated point

    !! Local variables
    real, dimension(3) :: v1
    integer :: i
    real    :: pi, scyl, sfac, tot, diddle

    pi = M_PI_sp

    scyl = 2.0*pi*r*al
    sfac = pi*r*r
    tot = scyl + sfac + sfac

    diddle = get_rand() * tot

    if (diddle < scyl) then
       call dovercylinder(c,n,r,al,t,p1)
    else if (diddle < scyl+sfac) then
       do i = 1,3
          v1(i) = c(i) + al*n(i)/2.0
       end do
       call doverdisk(v1,n,r,t,p1)
    else
       do i = 1,3
          v1(i) = c(i) - al*n(i)/2.0
       end do
       call doverdisk(v1,n,r,t,p1)
    end if

  end subroutine doversolcyl

  !! ----------------------------------------------------------------

  !> Generate a uniformly distributed point on a torus
  !! 
  !! @warning should be renamed
  !! 
  !! @warning uses anonymous numeric constants

  subroutine dovertorus(c, n, r1, r2, t, p1)

    !!GCC$ ATTRIBUTES unused :: n

    !! Declare arguments
    real, dimension(3)   :: c, n     !> Torus center and unit normal
    real, intent(in)     :: r1, r2   !> Torus radii
    real, dimension(3,3) :: t        !> Special rotation matrix
    real, dimension(3)   :: p1       !> Generated point

    !! Local variables
    real, dimension(3) ::p2
    real :: pi, phi, probkeep, theta, rmax, rstretch, rtest

    !! Get around compiler warning for unused dummy argument(s)
    real, dimension(3) :: n_2

    !! HACK !!
    if (.false.) then
       n_2(1) = n(1)
    end if

    pi = M_PI_sp

    do
       theta = get_rand()*2.0*pi

       !> Generate a point on a circle in the x-z plane
       p2(1) = r2*cos(theta) + r1
       p2(2) = 0.0
       p2(3) = r2*sin(theta)

       !> Record stretching data
       rstretch = p2(1)
       rmax = r1 + r2
       probkeep = rstretch/rmax

       !> Rotate about z-axis through random angle phi
       phi = get_rand()*2.0*pi

       p1(1) = cos(phi)*p2(1) - sin(phi)*p2(2)
       p1(2) = sin(phi)*p2(1) + cos(phi)*p2(2)
       p1(3) = p2(3)

       rtest = get_rand()

       if (rtest <= probkeep) exit
    end do

    call backtransform(c,t,p1)

  end subroutine dovertorus

  !! ----------------------------------------------------------------

  !> Distribute a point randomly over ellipsoid
  !! 
  !! @warning should be renamed
  !! 
  !! @warning spaghetti code!

  subroutine doverellipsoid(c, n1, n2, aa, bb, cc, t, p1)

    !!GCC$ ATTRIBUTES unused :: n1, n2

    use zeno_warnings, only : nell, rerr

    !! Declare arguments
    real, dimension(3)   :: c, n1, n2    !> Ellipsoid center and normals
    real, intent(in)     :: aa, bb, cc   !> Ellipsoid axial lengths
    real, dimension(3,3) :: t            !> Special rotation matrix
    real, dimension(3)   :: p1           !> Generated point

    !! Local variables
    real, dimension(3) :: p2, am, ak
    real :: as, bs, cs, phi, probkeep, rtest, theta
    real :: stretch, stretch1, stretch2, stretchmax

    !! Get around compiler warning for unused dummy argument(s)
    real, dimension(3) :: n1_2, n2_2

    !! HACK !!
    if (.false.) then
       n1_2(1) = n1(1)
       n2_2(2) = n2(2)
    end if

    !> Sort the eigenvalues:
    as = aa
    bs = bb
    cs = cc
    call sort3(as,bs,cs)    !  as > bs > cs
    stretchmax = (as*bs)**2 + ((as**2 - bs**2)**2)/4.0
    stretchmax = sqrt(stretchmax)

    do
       call sphere(p2,1.0)
       p1(1) = aa*p2(1)
       p1(2) = bb*p2(2)
       p1(3) = cc*p2(3)
       call makepolar(p2,theta,phi)
       am(1) = cos(theta)*cos(phi)
       am(2) = cos(theta)*sin(phi)
       am(3) = -sin(theta)
       ak(1) = -sin(phi)
       ak(2) = cos(phi)
       ak(3) = 0.0
       am(1) = aa*am(1)
       am(2) = bb*am(2)
       am(3) = cc*am(3)
       ak(1) = aa*ak(1)
       ak(2) = bb*ak(2)
       ak(3) = cc*ak(3)
       stretch1 = am(1)**2 + am(2)**2 + am(3)**2
       stretch2 = ak(1)**2 + ak(2)**2 + ak(3)**2
       stretch = sqrt(stretch1*stretch2)

       if (stretch > stretchmax) then
          if (nell == 0) then
             rerr = stretch/stretchmax
             nell = 1
          else if (nell == 1) then
             rerr = amax1(rerr,stretch/stretchmax)
          end if
       end if

       probkeep = stretch/stretchmax

       rtest = get_rand()
       if (rtest <= probkeep) exit
    end do

    call backtransform(c,t,p1)

  end subroutine doverellipsoid

  !! ----------------------------------------------------------------

  !> Generate a uniformly distributed point on a pillar.
  !! 
  !! @warning should be renamed

  subroutine doverpillar(x, y, zl, zh, p1)

    !! Declare arguments
    real, dimension(3) :: x, y, zl, zh    !> Pillar specs
    real, dimension(3) :: p1              !> Point generated

    !! local variables
    real, dimension(3) :: v1, v2, v3, v4
    real :: sa, bottom, top, side12, side13, side23
    real :: piddle

    call pillarsurf(x,y,zl,zh,sa,bottom,top,side12,side13,side23)

    piddle = get_rand()*sa

    if (piddle < bottom) then
       !> Find a point on the bottom triangle
       v1(1) = x(1)
       v1(2) = y(1)
       v1(3) = zl(1)
       v2(1) = x(2)
       v2(2) = y(2)
       v2(3) = zl(2)
       v3(1) = x(3)
       v3(2) = y(3)
       v3(3) = zl(3)
       call dovertriangle(v1,v2,v3,p1)

    else if (piddle < bottom+top) then
       !> Find a point on the top triangle
       v1(1) = x(1)
       v1(2) = y(1)
       v1(3) = zh(1)
       v2(1) = x(2)
       v2(2) = y(2)
       v2(3) = zh(2)
       v3(1) = x(3)
       v3(2) = y(3)
       v3(3) = zh(3)
       call dovertriangle(v1,v2,v3,p1)

    else if (piddle < bottom+top+side12) then
       !> Find a point on side12
       v1(1) = x(1)
       v1(2) = y(1)
       v1(3) = zl(1)
       v2(1) = x(2)
       v2(2) = y(2)
       v2(3) = zl(2)
       v3(1) = x(1)
       v3(2) = y(1)
       v3(3) = zh(1)
       v4(1) = x(2)
       v4(2) = y(2)
       v4(3) = zh(2)
       call dovertrap(v1,v2,v3,v4,p1)
          
    else if (piddle < bottom+top+side12+side13) then
       !> Find a point on side13
       v1(1) = x(1)
       v1(2) = y(1)
       v1(3) = zl(1)
       v2(1) = x(3)
       v2(2) = y(3)
       v2(3) = zl(3)
       v3(1) = x(1)
       v3(2) = y(1)
       v3(3) = zh(1)
       v4(1) = x(3)
       v4(2) = y(3)
       v4(3) = zh(3)
       call dovertrap(v1,v2,v3,v4,p1)

    else
       !> Find a point on side23
       v1(1) = x(2)
       v1(2) = y(2)
       v1(3) = zl(2)
       v2(1) = x(3)
       v2(2) = y(3)
       v2(3) = zl(3)
       v3(1) = x(2)
       v3(2) = y(2)
       v3(3) = zh(2)
       v4(1) = x(3)
       v4(2) = y(3)
       v4(3) = zh(3)
       call dovertrap(v1,v2,v3,v4,p1)

    end if

  end subroutine doverpillar

  !! ----------------------------------------------------------------

  !> Generate a uniformly distributed point on a trapezoid?
  !! 
  !! @warning should be renamed

  subroutine dovertrap(v1, v2, v3, v4, p1)

    !! Declare arguments
    real, dimension(3) :: v1, v2, v3, v4    !> Trapezoid vertices
    real, dimension(3) :: p1                !> Point generated

    !! Local variables
    real :: pearl, s13, s24

    s24 = v4(3) - v2(3)
    s13 = v3(3) - v1(3)

    pearl = get_rand()*(s24+s13)

    if (pearl < s13) then
       call dovertriangle(v1,v3,v4,p1)
    else
       call dovertriangle(v1,v2,v4,p1)
    end if

  end subroutine dovertrap

  !! ================================================================

end module zeno_shape_surface_point

!!! ================================================================

!!! Local Variables:
!!! mode: f90
!!! time-stamp-line-limit: 30
!!! fill-column: 72
!!! End:
