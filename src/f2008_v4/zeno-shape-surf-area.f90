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
!! This file provides the `zeno_shape_surface_area` module which defines
!! routines that compute the surface areas of various shapes.
!! 
!! @warning Need definitions for all shapes
!! 
!! @todo Code crying out for an object-oriented aproach with a `shape`
!! class hierarchy.
!!
!! @author Mark Mansfield (formerly of Stevens Institute) and students
!! @author Walid Keyrouz (walid.keyrouz@nist.gov)
!! @date   Fri Oct 25 09:33:26 2013 EDT
!! 
!! @defgroup zeno_shape_surface_area Shape Surface Areas
!! 
!! Groups routines that compute the surface areas of various shapes.
!! 
!! @warning should be merged as methods in a `shape` class hierarchy.
!
! Time-stamp: <2015-01-26 12:30:25 wtk>
! 
! ================================================================

module zeno_shape_surface_area

  !! ================================================================

  use numeric_constants

  use zeno_vectors, only : vector_difference, dotproduct, pythag0
  use zeno_sort_utils, only : sort3

  implicit none
  private

  public carea, pillarsurf

contains

  !! ================================================================
  !!
  !! Public routines
  !!
  !! ================================================================

  !> Tentative: computes surface area of shape objects.
  !!
  !! @warning may not handle all shape types!
  !!
  !! @warning uses anonymous numeric constants!

  subroutine carea(maxelts, eltype, bv, nelts, saar, total)

    use zeno_codes_data

    !! Declare arguments
    integer, intent(in)         :: maxelts    !> Array capacity
    integer, dimension(maxelts) :: eltype     !> Elements array
    real, dimension(maxelts,12) :: bv         !> TBD
    integer, intent(in)         :: nelts      !> Array size
    real, dimension(maxelts)    :: saar       !> TBD
    real, intent(out)           :: total      !> Total surface area?

    !! Local variables
    integer :: i
    real    :: sa

    total = 0.0

    do_all_elts: do i = 1,nelts

       call shape_surf(maxelts, eltype(i), i, bv, sa)

       saar(i) = sa
       total = total + sa

    end do do_all_elts

  end subroutine carea

  !! ----------------------------------------------------------------

  !> Compute the surface area of a pillar along with the areas of the
  !! individual faces.

  subroutine pillarsurf(x, y, zl, zh, sa, bottom, top, &
       &                side12, side13, side23)

    !! Declare arguments
    real, dimension(3) :: x, y, zl, zh    !> Pillar specs
    real, intent(out)  :: sa              !> Total surface area
    real, intent(out)  :: bottom, top     !> Area of bottom and top faces
    real, intent(out)  :: side12          !> Area of face `side12`
    real, intent(out)  :: side13          !> Area of face `side13`
    real, intent(out)  :: side23          !> Area of face `side23`

    !! Local variables
    real, dimension(3) :: v1, v2, v3
    real               :: alt, alt2, s1, s2, s3

    !> Use trisurf to get the areas of the top and bottom faces

    v1(1) = x(1)
    v1(2) = y(1)
    v1(3) = zl(1)
    v2(1) = x(2)
    v2(2) = y(2)
    v2(3) = zl(2)
    v3(1) = x(3)
    v3(2) = y(3)
    v3(3) = zl(3)
    call trisurf(v1,v2,v3,bottom)

    v1(1) = x(1)
    v1(2) = y(1)
    v1(3) = zh(1)
    v2(1) = x(2)
    v2(2) = y(2)
    v2(3) = zh(2)
    v3(1) = x(3)
    v3(2) = y(3)
    v3(3) = zh(3)
    call trisurf(v1,v2,v3,top)

    alt2 = (x(1)-x(2))**2 + (y(1)-y(2))**2
    alt = sqrt(alt2)
    s1 = zh(1) - zl(1)
    s2 = zh(2) - zl(2)
    side12 = alt*(0.5)*(s1+s2)

    alt2 = (x(1)-x(3))**2 + (y(1)-y(3))**2
    alt = sqrt(alt2)
    s1 = zh(1) - zl(1)
    s3 = zh(3) - zl(3)
    side13 = alt*(0.5)*(s1+s3)
    
    alt2 = (x(2)-x(3))**2 + (y(2)-y(3))**2
    alt = sqrt(alt2)
    s2 = zh(2) - zl(2)
    s3 = zh(3) - zl(3)
    side23 = alt*(0.5)*(s2+s3)

    sa = top + bottom + side12 + side13 + side23

  end subroutine pillarsurf

  !! ================================================================
  !!
  !! Private routines
  !!
  !! ================================================================

  !> Computes the surface area of one element.  Purpose of routine is to
  !! remove horrendeous type-related `if/switch` statement out of main
  !! loop.

  subroutine shape_surf(maxelts, elt_type, i, bv, sa)

    use zeno_codes_data

    !! Declare arguments
    integer, intent(in)         :: maxelts    !> Array capacity
    integer, intent(in)         :: elt_type   !> Element shape's type
    integer, intent(in)         :: i          !> Elt index
    real, dimension(maxelts,12) :: bv         !> TBD
    real, intent(out)           :: sa         !> Element shape's surf. area

    !! Local variables
    real, dimension(3) :: v1, v2, v3
    real, dimension(3) :: x, y, zl, zh
    integer :: j
    real    :: aa, al, bb, bottom, cc, pi, r, r1, r2
    real    :: side, side12, side13, side23, top

    pi = M_PI_sp

    Cs_Elt_Type: select case (elt_type)

    case (cube_code)
       side = bv(i,4)
       sa = 6.0*side*side

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
       call pillarsurf(x,y,zl,zh,sa,bottom,top,side12,side13,side23)

    case (sphere_code)
       r = bv(i,4)
       sa = 4.0 * pi * r * r

    case (triangle_code)
       do j = 1,3
          v1(j) = bv(i,j)
          v2(j) = bv(i,j+3)
          v3(j) = bv(i,j+6)
       end do
       call trisurf(v1,v2,v3,sa)

    case (disk_code)
       r = bv(i,7)
       sa = pi * r * r

    case (open_cylinder_code)
       r = bv(i,7)
       al = bv(i,8)
       sa = 2.0 * pi * al * r

    case (solid_cylinder_code)
       r = bv(i,7)
       al = bv(i,8)
       sa = 2.0 * pi * al * r
       sa = sa + 2.0 * pi * r * r

    case (donut_code)
       r1 = bv(i,7)
       r2 = bv(i,8)
       sa = 4.0 * pi * pi * r1 * r2

    case (ellipsoid_code)
       aa = bv(i,10)
       bb = bv(i,11)
       cc = bv(i,12)
       call ellsurf(aa,bb,cc,sa)

    end select Cs_Elt_Type


  end subroutine shape_surf

  !! ----------------------------------------------------------------

  !> Compute the surface area of a triangle

  subroutine trisurf(v1, v2, v3, sa)

    !! Declare arguments
    real, dimension(3) :: v1, v2, v3    !> Triangle vertices
    real, intent(out)  :: sa            !> Triangle area

    !! Local variables
    real, dimension(3) :: p23, p21, q, av
    integer :: i
    real    :: a, bot, h, top, tt

    call vector_difference(v2,v3,p23)
    call vector_difference(v2,v1,p21)
    call dotproduct(p23,p21,top) 
    call dotproduct(p21,p21,bot)
    tt = top/bot
    h = sqrt(bot)
    do i = 1,3
       q(i) = tt*v1(i) + (1.0-tt)*v2(i)
    end do
    call vector_difference(v3,q,av)
    call pythag0(av,a)
    sa = a*h/2.0

  end subroutine trisurf

  !! ----------------------------------------------------------------

  !> Compute surface area of ellipsoid
  !!
  !! @warning Uses anonymous numeric constants.

  subroutine ellsurf(aa, bb, cc, sa)

    !! Declare arguments
    real, intent(in)  :: aa, bb, cc    !> Ellipsoid axial lengths
    real, intent(out) :: sa            !> Ellipsoid surface area

    !! Local variables

    integer :: iphi, itheta
    real    :: as, bs, cs, hphi, htheta
    real    :: p, pi, phi, rr, sum, t1, t2, theta

    !> Sort the eigenvalues:
    as = aa
    bs = bb
    cs = cc
    call sort3(as,bs,cs)    !  as > bs > cs

    !> Evaluate surface integral by double simpsons

    pi = M_PI_sp

    !> @warning compares two pairs of reals with ==
    if (as == bs .and. bs == cs) then
       sa = 4.0 * pi * cs * cs
       return

       !> @warning compares two reals with ==
    else if (as == bs) then
       phi = asin (sqrt (1.0 - (cs/as)**2) )   
       t1 = alog(tan(pi/4.0 + phi/2.0))
       t1 = t1*as*cs*cs/sqrt(as**2 - cs**2)
       t2 = sqrt(as**2 - 2.0*cs**2 + cs**4/as**2)
       t2 = t2*as
       sa = 2.0*pi*(cs**2 + t1 + t2)
       return

       !> @warning compares two reals with ==
    else if (bs == cs) then
       phi = asin (sqrt (1.0 - (cs/as)**2) )   
       t1 = cs*as*as/sqrt(as**2 - cs**2)
       sa = 2.0*pi*(cs**2 + t1*phi)
       return
    end if

    sa = 0.0
    hphi = 2.0*pi/1000.0
    do iphi = 0,1000
       phi = float(iphi)*hphi
       sum = 0.0
       htheta = pi/1000.0
       do itheta = 0,1000
          theta = float(itheta)*htheta
          p = (bs*cs*sin(theta)*cos(phi))**2
          p = p + (as*cs*sin(theta)*sin(phi))**2
          p = p + (as*bs*cos(theta))**2
          p = sin(theta)*sqrt(p)
          if (itheta == 0) then
             rr = 1.0
          else if (itheta == 1000) then
             rr = 1.0
          else if (mod(itheta,2) == 1) then
             rr = 4.0
          else
             rr = 2.0
          end if
          sum = sum + rr*p
       end do
       sum = htheta*sum/3.0
       if (iphi == 0) then
          rr = 1.0
       else if (iphi == 1000) then
          rr = 1.0
       else if (mod(iphi,2) == 1) then
          rr = 4.0
       else
          rr = 2.0
       end if
       sa = sa + rr*sum
    end do
    sa = hphi*sa/3.0

  end subroutine ellsurf

  !! ================================================================

end module zeno_shape_surface_area

!!! ================================================================

!!! Local Variables:
!!! mode: f90
!!! time-stamp-line-limit: 30
!!! fill-column: 72
!!! End:
