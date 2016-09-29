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
!! This file provides the `zeno_shape_distances` module which defines
!! routines for computing distances to geometric shapes.
!! 
!! @warning Need definitions for all shapes
!! 
!! @todo Code crying out for an object-oriented aproach with a `shape`
!! class hierarchy.
!!
!! @author Mark Mansfield (formerly of Stevens Institute) and students
!! @author Walid Keyrouz (walid.keyrouz@nist.gov)
!! @date   Tue Oct 22 17:52:30 2013 EDT
!! 
!! @defgroup zeno_shape_distances Shape Distances
!! 
!! Groups together routines that compute point distances from the origin
!! to various shapes (minimum, maximum, etc.).
!!
!! @warning Need definitions for all shapes.
!! 
!! @todo Code crying out for an object-oriented aproach with a `shape`
!! class hierarchy.
! 
! Time-stamp: <2015-01-25 14:14:47 wtk>
! 
! ================================================================

module zeno_shape_distances

  !! ================================================================

  use zeno_vectors, only : pythag0, pythag, &
       & dotproduct, cross_product, scalar_product, vector_difference, &
       & rotate, backtransform
  use zeno_ellipsoid_utils, only : exellipsoid
  use zeno_sort_utils, only : listersort
  use zeno_sphere, only : sphere

  implicit none
  private 

  public &
       & distance, bridge, &
       & maxcube, procube, mincube, &
       & maxpillar, propillar, minpillar, &
       & maxsphere, prosphere, minsphere, &
       & maxtriangle, protriangle, mintriangle, &
       & maxdisk, prodisk, mindisk, &
       & maxcylinder, procylinder, mincylinder, &
       & maxtorus, protorus, mintorus, &
       & maxsolcyl, prosolcyl, minsolcyl, &
       & maxellipsoid, proellipsoid, minellipsoid, &
       & minedge

contains

  !! ================================================================
  !!
  !! Public routines
  !!
  !! ================================================================

  !> @brief Computes the distance from point `p` to the body's surface.
  !!
  !! Compute distance from point `p` to the surface of the body
  !! The distance is returned as `ds`.  The sphere centered at point
  !! `p` having radius `ds` is called the jump sphere.
  !! 
  !! @todo
  !! Document current version of routine which is now calls different
  !! routines depending on whether `bubble` is properly defined.
  !!
  !! @warning may be renamed to something like `dist_to_point` and
  !! turned into a method.
  !!
  !! @warning routine seems to be the "abstract" method gateway that
  !! invokes a bunch of other @c "min*" routines.  This may lead to
  !! having these routines becoming @e private ones.
  !! 
  !! @warning too many arguments

  subroutine distance(maxelts, eltype, bv, nelts, rotations, p, ds, &
       &              bubble, bubble_rad, nebtab, nneb, ninn, rlist, &
       &              nearto)

    use numeric_kinds
    use zeno_debug

    !! Declare arguments
    integer, intent(in)          :: maxelts    !> Elements array capacity
    integer, dimension(maxelts)  :: eltype     !> Elements array
    real, dimension(maxelts,12)  :: bv         !> TBD
    integer, intent(in)          :: nelts      !> Elements array size
    real, dimension(maxelts,3,3) :: rotations  !> Rotations array
    real, dimension(3)           :: p          !> Arbitrary point
    real, intent(out)            :: ds         !> Distance computed
    real, dimension(3)           :: bubble     !> Center of neighbors bubble
    real                         :: bubble_rad !> Radius of neighbors bubble
    integer, dimension(maxelts)  :: nebtab     !> Neighbors table
    integer                      :: nneb       !> Number of neighbors
    integer                      :: ninn       !> Number of elements
    real, dimension(maxelts)     :: rlist      !> Radius of launch sphere
    integer                      :: nearto     !> TBD

    !! Local variables
    real    :: close

    !! ----------------------------------------------------------------

    !> Update `log` data

    calls_DIST_cnt = calls_DIST_cnt + 1

    !> Compute distance from point `p` to the surface of the body.
    !! The distance is returned as `ds`.  The sphere centered at
    !! point `p` having radius `ds` is called the @e jump sphere.

    !! ----------------------------------------------------------------

    If_Bubble_Rad_ge_0: if (bubble_rad >= 0.0) then

       call distance_scan_neighbors(maxelts, eltype, bv, nelts, rotations, &
            &                       p, ds, nebtab, nneb, ninn, rlist, &
            &                       nearto, close, bubble, bubble_rad)

       !> @warning compares two reals with ==
       if (ds /= close) then
          print *,'ds,close = ',ds,close
          print *,' ds /= close    error'
          stop '2'
       end if

    else                        !! (bubble_rad < 0)

       call distance_scan_all_elts(maxelts, eltype, bv, nelts, rotations, &
            &                      p, ds, bubble, bubble_rad, nebtab, &
            &                      nneb, ninn, rlist, nearto, close)

       !> @warning compares two reals with ==
       if (ds /= close) then
          print *,'ds,close = ',ds,close
          print *,' ds /= close    error'
          stop '1'
       end if

    end if If_Bubble_Rad_ge_0

  end subroutine distance

  !! ----------------------------------------------------------------

  !> Scan elements in the `neighbors` table for the smallest distance.

  subroutine distance_scan_neighbors(maxelts, eltype, bv, nelts, rotations, &
       &                             p, ds, nebtab, nneb, ninn, rlist, &
       &                             nearto, close, bubble, bubble_rad)

    use numeric_kinds
    use zeno_debug

    !! Declare arguments
    integer, intent(in)          :: maxelts    !> Elements array capacity
    integer, dimension(maxelts)  :: eltype     !> Elements array
    real, dimension(maxelts,12)  :: bv         !> TBD
    integer, intent(in)          :: nelts      !> Elements array size
    real, dimension(maxelts,3,3) :: rotations  !> Rotations array
    real, dimension(3)           :: p          !> Arbitrary point
    real, intent(out)            :: ds         !> Distance computed
    integer, dimension(maxelts)  :: nebtab     !> Neighbors list
    integer                      :: nneb       !> Number of neighbors
    integer                      :: ninn       !> Actual num. of neighbors
    real, dimension(maxelts)     :: rlist      !> Radius of launch sphere
    integer                      :: nearto     !> TBD
    real, intent(in out)         :: close      !> TBD
    real, dimension(3)           :: bubble     !> Center of neighbors bubble
    real                         :: bubble_rad !> Radius of neighbors bubble

    !! Local variables
    integer :: i, ii
    real    :: d, pb

    !! ----------------------------------------------------------------
    
    !> Update `log` data
    calls_DIST_scan_neighbors = calls_DIST_scan_neighbors + 1

    !> The current bubble and neighborlist are assumed to be OK -- we
    !! will scan all elements if we later find out they are not

    !> Scan neighbor table for smallest distance

    !> Case 2: bubble_rad >= 0
    !! Generating the jump sphere by probing nearby body elements.  If
    !! any part of the jump sphere is later found to lie outside the
    !! bubble, we will redo this step.

    nearto = 0

    Do_Ninn : do ii = 1,ninn

       i = nebtab(ii)  

       call min_shape(eltype(i), i, maxelts, bv, rotations, p, d)

       if (ii == 1) then
          ds = d
       else
          ds = amin1(ds,d)
       end if

       if (nearto == 0) then
          close = d
          nearto = i
       else
          if (d < close) then
             close = d
             nearto = i
          end if
       end if

    end do Do_Ninn
    
    !! ----------------------------------------------------------------
    
    !> At this point we have identified a jump sphere of radius `ds`
    !! centered on the point `p`, which is the current position of the
    !! random walker.  However, it was determined only by probing
    !! the nearby body elements.  Therefore, it must be rejected
    !! if any part of this sphere lies outside of the bubble generated
    !! at the time that this particular list of nearby elements
    !! was generated.  So if it is, we jump back up and generate the
    !! jump sphere by probing all body elements.
    
    !> `pb` is the distance from point p to the center of the bubble
    pb = 0.0
    do i = 1,3
       pb = pb + (bubble(i)-p(i))**2
    end do
    pb = sqrt(pb)

    if (pb+ds > bubble_rad) then

       !> Update `log` data
       calls_DIST_re_comp_bbl = calls_DIST_re_comp_bbl + 1

       call distance_scan_all_elts(maxelts, eltype, bv, nelts, rotations, &
            &                      p, ds, bubble, bubble_rad, nebtab, &
            &                      nneb, ninn, rlist, nearto, close)

       !> @warning compares two reals with ==
       if (ds /= close) then
          print *,'ds,close = ',ds,close
          print *,' ds /= close    error'
          stop '1'
       end if

    end if

  end subroutine distance_scan_neighbors

  !! ----------------------------------------------------------------

  !> @brief Computes min distance by scanning all elements.
  !! 
  !! Compute min distance by scanning all elements.  Also, computes the
  !! bubble's radius and number of neighbors in it.  Purpose of routine
  !! is to remove `goto` statements and avoid repeating code in
  !! subroutine `distance`!
  !! 
  !! @warning argument list is too long!
  !! 
  !! @todo reimplement using a spatial data structure.

  subroutine distance_scan_all_elts(maxelts, eltype, bv, nelts, rotations, &
       &                            p, ds, bubble, bubble_rad, nebtab, &
       &                            nneb, ninn, rlist, nearto, close)

    use numeric_kinds
    use zeno_debug

    !! Declare arguments
    integer, intent(in)          :: maxelts    !> Elements array capacity
    integer, dimension(maxelts)  :: eltype     !> Elements array
    real, dimension(maxelts,12)  :: bv         !> TBD
    integer, intent(in)          :: nelts      !> Elements array size
    real, dimension(maxelts,3,3) :: rotations  !> Rotations array
    real, dimension(3)           :: p          !> Arbitrary point
    real, intent(out)            :: ds         !> Distance computed
    real, dimension(3)           :: bubble     !> Bubble center?
    real                         :: bubble_rad !> Bubble radius?
    integer, dimension(maxelts)  :: nebtab     !> Neighbors list
    integer, intent(in)          :: nneb       !> Max num. of neighbors
    integer, intent(in out)      :: ninn       !> Actual num. of neighbors
    real, dimension(maxelts)     :: rlist      !> TBD
    integer                      :: nearto     !> TBD
    real                         :: close      !> TBD

    !! Local variables
    integer :: i
    real    :: d

    !! ----------------------------------------------------------------

    !> Update `log` data
    calls_DIST_scan_all_elts = calls_DIST_scan_all_elts + 1

    !! ----------------------------------------------------------------

    ninn = 0
    nearto = 0
    
    Do_Nelts : do i = 1,nelts

       call min_shape(eltype(i), i, maxelts, bv, rotations, p, d)

       !> `listersort` is called after every `min_shape` invocation
       !> to keep at most `nneb` neighbors in sorted list.

       ninn = ninn + 1
       rlist(ninn) = d
       nebtab(ninn) = i
       call listersort(rlist,ninn,nneb,nebtab,maxelts)

       if (nearto == 0) then
          nearto = i
          close = d
       else
          if (d < close) then
             close = d
             nearto = i
          end if
       end if

    end do Do_Nelts

    ds = rlist(1)
    bubble_rad = rlist(ninn)
    do i = 1,3
       bubble(i) = p(i)
    end do

  end subroutine distance_scan_all_elts

  !! ----------------------------------------------------------------

  !> @brief Tentative: computes projection of all shapes onto a
  !! particular axis.
  !! 
  !! Computes projection of all shapes onto a particular axis.
  !! Upon entry, `v` is an arbitrary unit vector.  `d` is returned
  !! as the span of the body projected onto the vector `v`.
  !!
  !! @warning Implementation looks suspicious because it only invokes @c
  !! pro_shape.

  subroutine bridge(maxelts, eltype, bv, nelts, rotations, v, d)

    use zeno_codes_data

    !! Declare arguments
    integer, intent(in)          :: maxelts   !> Array capacity
    integer, dimension(maxelts)  :: eltype    !> Elements array
    real, dimension(maxelts,12)  :: bv        !> Elt. body values
    integer, intent(in)          :: nelts     !> Array size
    real, dimension(maxelts,3,3) :: rotations !> Array of elt. rotations
    real, dimension(3)           :: v         !> Projection vector
    real, intent(out)            :: d         !> Projection span

    !! Local variables
    integer :: i
    real    :: d1, d2, s1, s2

    Do_Nelts : do i = 1,nelts

       call pro_shape(eltype(i), i, maxelts, bv, rotations, v, s1, s2)

       if (i == 1) then
          d1 = s1
          d2 = s2
       else
          d1 = amin1(d1,s1)
          d2 = amax1(d2,s2)
       end if

    end do Do_Nelts

    d = d2 - d1
        
  end subroutine bridge

  !! ================================================================

  !> Given a cube with low corner at `v1` and high corner at `v2`,
  !! return in `dd`, the maximum distance between all points on the
  !! surface of the cube and the origin.

  subroutine maxcube(v1, side, dd)

    !! Declare arguments
    real, dimension(3) :: v1
    real, intent(in)   :: side
    real, intent(out)  :: dd

    !! Local variables

    integer :: ix, iy, iz
    real    :: d
    real, dimension(3) :: v2

    do ix = 0,1
       do iy = 0,1
          do iz = 0,1
             v2(1) = v1(1) + float(ix)*side
             v2(2) = v1(2) + float(iy)*side
             v2(3) = v1(3) + float(iz)*side
             call pythag0(v2,d)
             if (ix+iy+iz == 0) then
                dd = d
             else
                dd = amax1(dd,d)
             end if
          end do
       end do
    end do

  end subroutine maxcube

  !! ----------------------------------------------------------------

  !> Projects a shape element onto a particular vector.  Purpose of
  !> routine is to move horrendous `if/switch` statement on `eltype`
  !> out of loop over all or a subset of elements
  !>
  !> @warning rename to something more meaningful.

  subroutine pro_shape(elt_type, i, maxelts, bv, rotations, v, s1, s2)

    use zeno_codes_data

    !! Declare arguments
    integer, intent(in)          :: elt_type  !> Element type
    integer, intent(in)          :: i         !> Element index
    integer, intent(in)          :: maxelts   !> Array capacity
    real, dimension(maxelts,12)  :: bv        !> Elt. body values
    real, dimension(maxelts,3,3) :: rotations !> Array of elt. rotations
    real, dimension(3)           :: v         !> Projection vector
    real, intent(out)            :: s1, s2    !> Projection span

    !! Local arguments
    real, dimension(3)   :: c, v1, v2, v3, n, n1, n2
    real,  dimension(3)  :: x, y, zl, zh
    real, dimension(3,3) :: t
    integer :: j, k
    real    :: aa, al, bb, cc, r, r1, r2, side
    
    Case_Elt_Type_i: select case (elt_type)

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
       call propillar(x,y,zl,zh,v,s1,s2)

    case(cube_code)
       do j = 1,3
          v1(j) = bv(i,j)
       end do
       side = bv(i,4)
       call procube(v1,side,v,s1,s2)

    case (sphere_code)
       do j = 1,3
          c(j) = bv(i,j)
       end do
       r = bv(i,4)
       call prosphere(c,r,v,s1,s2)

    case(triangle_code)
       do j = 1,3
          v1(j) = bv(i,j)
          v2(j) = bv(i,j+3)
          v3(j) = bv(i,j+6)
       end do
       call protriangle(v1,v2,v3,v,s1,s2)

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
       call prodisk(c,n,r,t,v,s1,s2)

    case(open_cylinder_code)
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
       call procylinder(c,n,r,al,t,v,s1,s2)

    case(solid_cylinder_code)
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
       call prosolcyl(c,n,r,al,t,v,s1,s2)

    case(donut_code)
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
       call protorus(c,n,r1,r2,t,v,s1,s2)

    case(ellipsoid_code)
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
       call proellipsoid(c,n1,n2,aa,bb,cc,t,v,s1,s2)

    end select Case_Elt_Type_i

  end subroutine pro_shape

  !! ----------------------------------------------------------------

  !> @b Original-comment: Given a cube with lower corner at `v1`, 
  !! `side`, and an arbitrary direction `v`, determine the distance
  !! between the two tangents (what tangents?).
  !! 
  !! @warning what tangents and what is a procube?
  
  subroutine procube(v1, side, v, s1, s2)

    !! Declare arguments
    real, dimension(3) :: v1        !> TBD
    real, intent(in)   :: side      !> TBD
    real, dimension(3) :: v         !> TBD
    real               :: s1        !> TBD
    real               :: s2        !> TBD

    !! Local variables
    integer :: ix, iy, iz, it
    real    ::    dot
    real, dimension(3) :: vt

    it = 0

    do ix = 0,1
       do iy = 0,1
          do iz = 0,1
             vt(1) = v1(1) + float(ix)*side
             vt(2) = v1(2) + float(iy)*side
             vt(3) = v1(3) + float(iz)*side
             call dotproduct(v,vt,dot)
             if (it == 0) then
                s1 = dot
                s2 = dot
                it = 1
             else
                s1 = amin1(s1,dot)
                s2 = amax1(s2,dot)
             end if
          end do
       end do
    end do

  end subroutine procube

  !! ----------------------------------------------------------------

  !> Computes minimum distance of a point to a particular shape.
  !! Purpose of routine is to move horrendous `if/switch` statement on
  !! `eltype(i)` out of loop over all or a subset of elements.
  !!
  !! @warning rename to something more meaningful.

  subroutine min_shape(elt_type, i, maxelts, bv, rotations, p, d)

    use zeno_codes_data

    !! Declare arguments
    integer, intent(in)          :: elt_type   !> Element type
    integer, intent(in)          :: i          !> Element index
    integer, intent(in)          :: maxelts    !> Elements array capacity
    real, dimension(maxelts,12)  :: bv         !> TBD
    real, dimension(maxelts,3,3) :: rotations  !> Rotations array
    real, dimension(3)           :: p          !> Arbitrary point
    real, intent(out)            :: d

    !! Local variables
    real, dimension(3) :: c, v1, v2, v3, n, n1, n2
    real, dimension(3) :: x, y, zl, zh
    real, dimension(3,3) :: t
    integer :: j, k
    real :: aa, al, bb, cc
    real :: r, r1, r2, side

    Cs_Elt_Type_Nelts : select case (elt_type)

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
       call minpillar(x,y,zl,zh,p,d)

    case (cube_code)
       do j = 1,3
          v1(j) = bv(i,j)
       end do
       side = bv(i,4)
       call mincube(v1,side,p,d)

    case (sphere_code)
       do j = 1,3
          c(j) = bv(i,j)
       end do
       r = bv(i,4)
       call minsphere(c,r,p,d)

    case (triangle_code)
       do j = 1,3
          v1(j) = bv(i,j)
          v2(j) = bv(i,j+3)
          v3(j) = bv(i,j+6)
       end do
       call mintriangle(v1,v2,v3,p,d)

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
       call mindisk(c,n,r,t,p,d)

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
       call mincylinder(c,n,r,al,t,p,d)

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
       call minsolcyl(c,n,r,al,t,p,d)

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
       call mintorus(c,n,r1,r2,t,p,d)

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
       call minellipsoid(c,n1,n2,aa,bb,cc,t,p,d)

    end select Cs_Elt_Type_Nelts

  end subroutine min_shape

  !! ----------------------------------------------------------------

  !> Given a cube with its low corner at `v1`, with a side of 
  !! `side`, and an arbitrary point `p` outside the cube, return the
  !! minimum distance, `d`, from `p` to the cube.

  subroutine mincube(v1, side, p, d)

    use zeno_warnings, only : ncube, ferr

    !! Declare arguments
    real, dimension(3) :: v1        !> cube lower corner
    real, intent(in)   :: side      !> cube side length
    real, dimension(3) :: p         !> point outside the cube
    real, intent(out)  :: d         !> min distance from p to cube

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

    if (pr(3) < s2) then
       if (ncube == 0) then
          ncube = 1
          ferr = pr(3)/s2
       else
          ferr = amin1(ferr,pr(3)/s2)
       end if
       d = 0.0
       return
    end if

    if (pr(2) < s2) then
       d = pr(3) - s2
       return
    end if

    if (pr(1) < s2) then
       d = (pr(3)-s2)**2 + (pr(2)-s2)**2
       d = sqrt(d)
       return
    end if

    d = (pr(1)-s2)**2 + (pr(2)-s2)**2 + (pr(3)-s2)**2
    d = sqrt(d)

  end subroutine mincube

  !! ----------------------------------------------------------------

  !> Maximum distance from the origin to any point on the surface of a
  !! pillar.
  !!
  !! Because pillars are convex, you just need to check six vertices.
  !!
  !! @warning what is a pillar?
 
  subroutine maxpillar(x, y, zl, zh, dd)

    !! Declare arguments
    real, dimension(3) :: x    !> TBD
    real, dimension(3) :: y    !> TBD
    real, dimension(3) :: zl   !> TBD
    real, dimension(3) :: zh   !> TBD
    real, intent(out)  :: dd   !> Maximum distance

    !! Local variables
    integer :: i, is
    real    :: d
    real, dimension(3) :: v

    is = 0

    do i = 1,3
       v(1) = x(i)
       v(2) = y(i)
       v(3) = zl(i)
       call pythag0(v,d)
       if (is == 0) then
          is = 1
          dd = d
       else
          dd = amax1(dd,d)
       end if
       v(3) = zh(i)
       call pythag0(v,d)
       dd = amax1(dd,d)
    end do

  end subroutine maxpillar

  !! ----------------------------------------------------------------
  
  !> Project a pillar surface onto a line.
  !! 
  !! @warning what is a pillar?
  
  subroutine propillar(x, y, zl, zh, v, s1, s2)

    !! Declare arguments
    real, dimension(3) :: x       !> TBD
    real, dimension(3) :: y       !> TBD
    real, dimension(3) :: zl      !> TBD
    real, dimension(3) :: zh      !> TBD
    real, dimension(3) :: v       !> Direction vector of line
    real, intent(out)  :: s1      !> Projection span, point 1
    real, intent(out)  :: s2      !> Projection span, point 2

    !! Local variables

    integer :: i, is
    real    :: dot
    real, dimension(3) :: p

    is = 0

    do i = 1,3
       p(1) = x(i)
       p(2) = y(i)
       p(3) = zl(i)
       call dotproduct(v, p, dot)
       if (is == 0) then
          s1 = dot
          s2 = dot
          is = 1
       else
          s1 = amin1(s1, dot)
          s2 = amax1(s2, dot)
       end if
       p(3) = zh(i)
       call dotproduct(v,p,dot)
       s1 = amin1(s1,dot)
       s2 = amax1(s2,dot)
    end do

  end subroutine propillar

  !! ----------------------------------------------------------------

  !> Minimum distance from point `p` to any point on a pillar suface.
  !!  
  !! @warning what is a minpillar?

  subroutine minpillar(x, y, zl, zh, p, d)

    !! Declare arguments
    real, dimension(3) :: x         !> TBD
    real, dimension(3) :: y         !> TBD
    real, dimension(3) :: zl        !> TBD
    real, dimension(3) :: zh        !> TBD
    real, dimension(3) :: p         !> Arbitrary point
    real, intent(out)  :: d         !> Minimum distance

    !! Local variables
    real :: d1
    real, dimension(3) :: v1, v2, v3

    !> First triangle will be the bottom face

    v1(1) = x(1)
    v1(2) = y(1)
    v1(3) = zl(1)
    v2(1) = x(2)
    v2(2) = y(2)
    v2(3) = zl(2)
    v3(1) = x(3)
    v3(2) = y(3)
    v3(3) = zl(3)

    call mintriangle(v1,v2,v3,p,d1)
    d = d1

    !> Second triangle will be the top face
    
    v1(1) = x(1)
    v1(2) = y(1)
    v1(3) = zh(1)
    v2(1) = x(2)
    v2(2) = y(2)
    v2(3) = zh(2)
    v3(1) = x(3)
    v3(2) = y(3)
    v3(3) = zh(3)
    call mintriangle(v1,v2,v3,p,d1)
    d = amin1(d1,d)

    v1(1) = x(1)
    v1(2) = y(1)
    v1(3) = zl(1)
    v2(1) = x(2)
    v2(2) = y(2)
    v2(3) = zl(2)
    v3(1) = x(1)
    v3(2) = y(1)
    v3(3) = zh(1)
    call mintriangle(v1,v2,v3,p,d1)
    d = amin1(d1,d)

    v1(1) = x(1)
    v1(2) = y(1)
    v1(3) = zh(1)
    v2(1) = x(2)
    v2(2) = y(2)
    v2(3) = zh(2)
    v3(1) = x(2)
    v3(2) = y(2)
    v3(3) = zl(2)
    call mintriangle(v1,v2,v3,p,d1)
    d = amin1(d1,d)

    v1(1) = x(1)
    v1(2) = y(1)
    v1(3) = zl(1)
    v2(1) = x(3)
    v2(2) = y(3)
    v2(3) = zl(3)
    v3(1) = x(1)
    v3(2) = y(1)
    v3(3) = zh(1)
    call mintriangle(v1,v2,v3,p,d1)
    d = amin1(d1,d)

    v1(1) = x(1)
    v1(2) = y(1)
    v1(3) = zh(1)
    v2(1) = x(3)
    v2(2) = y(3)
    v2(3) = zh(3)
    v3(1) = x(3)
    v3(2) = y(3)
    v3(3) = zl(3)
    call mintriangle(v1,v2,v3,p,d1)
    d = amin1(d1,d)

    v1(1) = x(2)
    v1(2) = y(2)
    v1(3) = zl(2)
    v2(1) = x(3)
    v2(2) = y(3)
    v2(3) = zl(3)
    v3(1) = x(2)
    v3(2) = y(2)
    v3(3) = zh(2)
    call mintriangle(v1,v2,v3,p,d1)
    d = amin1(d1,d)

    v1(1) = x(2)
    v1(2) = y(2)
    v1(3) = zh(2)
    v2(1) = x(3)
    v2(2) = y(3)
    v2(3) = zh(3)
    v3(1) = x(3)
    v3(2) = y(3)
    v3(3) = zl(3)
    call mintriangle(v1,v2,v3,p,d1)
    d = amin1(d1,d)

  end subroutine minpillar

  !! ----------------------------------------------------------------

  !> Given a sphere with center `c` and radius `r`, return in `dd`
  !! the maximum distance between all points on the surface of the
  !! sphere and the origin.

  subroutine maxsphere(c, r, dd)

    !! Declare arguments
    real, dimension(3) :: c
    real, intent(in)   :: r
    real, intent(out)  :: dd

    !! Local variables
    real :: d1

    call pythag0(c,d1)
    dd = d1 + r

  end subroutine maxsphere

  !! ----------------------------------------------------------------

  !> Given a sphere with center at `c` and radius `r`, and a projection
  !! axis `v`, return in `s1` and `s2` the span of the projection (`s1`
  !! is the beginning and `s2` is the end).
  !!
  !! @warning assumptions re point being inside or outside sphere?

  subroutine prosphere(c, r, v, s1, s2)

    !! Declare arguments
    real, dimension(3) :: c         !> Center of sphere?
    real, intent(in)   :: r         !> Radius of sphere
    real, dimension(3) :: v         !> Projection axis
    real, intent(out)  :: s1        !> Projection span start
    real, intent(out)  :: s2        !> Projection span finish

    !! Local variables
    real :: dot

    call dotproduct(c, v, dot)
    s1 = dot - r
    s2 = dot + r

  end subroutine prosphere

  !! ----------------------------------------------------------------

  !> Given a sphere with center at `c` and radius `r`, and an
  !! arbitrary point `p`, return in `d` the minimum distance between
  !! all points on the sphere and the point `p`.
  !!
  !! @warning assumptions re point being inside or outside sphere?

  subroutine minsphere(c, r, p, d)

    !! Declare arguments
    real, dimension(3) :: c         !> Center of sphere
    real, intent(in)   :: r         !> Radius of sphere
    real, dimension(3) :: p         !> Arbitrary point
    real, intent(out)  :: d         !> Minimum distance

    !! Local variables
    real :: r1

    call pythag(c,p,r1)
    d = abs(r1-r)

  end subroutine minsphere

  !! ----------------------------------------------------------------

  !> Given a triangle with vertices (`v1`, `v2`, `v3`), return the
  !! maximum distance between any point in the triangle to the origin
  !! in `dd`.
  !! 
  !! @warning Assumes one of the three vertices is at the maximum.
 
  subroutine maxtriangle(v1, v2, v3, dd)

    !! Declare arguments
    real, dimension(3) :: v1        !> Triangle vertex v1
    real, dimension(3) :: v2        !> Triangle vertex v2
    real, dimension(3) :: v3        !> Triangle vertex v3
    real, intent(out)  :: dd        !> Maxmum distance

    !! Local variables
    real :: d1, d2, d3, dd1

    call pythag0(v1,d1)
    call pythag0(v2,d2)
    call pythag0(v3,d3)
    dd1 = amax1(d1,d2)
    dd = amax1(dd1,d3)

  end subroutine maxtriangle

  !! ----------------------------------------------------------------

  !> @brief Project a triangle onto an arbitary line
  !!
  !! Given a triangle with vertices (`v1`, `v2`, `v3`) and an
  !! arbitrary line `v`, return in `s1` and `s2` the projection
  !! span of the triangle onto the line (`origin`, `v`).

  subroutine protriangle(v1, v2, v3, v, s1, s2)

    !! Declare arguments
    real, dimension(3) :: v1, v2, v3     !> Triangle vertices
    real, dimension(3) :: v              !> Direction vector of line
    real, intent(out)  :: s1, s2         !> Projection span

    !! Local variables
    real :: dot

    call dotproduct(v1,v,dot)
    s1 = dot
    s2 = dot
    call dotproduct(v2,v,dot)
    s1 = amin1(s1,dot)
    s2 = amax1(s2,dot)
    call dotproduct(v3,v,dot)
    s1 = amin1(s1,dot)
    s2 = amax1(s2,dot)

  end subroutine protriangle

  !! ----------------------------------------------------------------

  !> Given a triangle with vertices (`v1`, `v2`, `v3`) and an
  !! arbitrary point `p`, return in `d` the minimum distance between
  !! all points in the triangle and the point `p`.

  subroutine mintriangle(v1, v2, v3, p, d)

    !! Declare arguments
    real, dimension(3) :: v1, v2, v3     !> Triangle vertices
    real, dimension(3) :: p              !> Arbitrary point
    real, intent(out)  :: d              !> Minimum distance

    !! Local variables
    real, dimension(3) :: b1, b2, b3
    real, dimension(3) :: p1, t3b3, q
    real :: d1, d11, d12, d2, d22, d3, db3b3, dp1b3, dx
    real :: q1, q2, t1, t2, t3

    call vector_difference(v2,v1,b1)
    call vector_difference(v3,v1,b2)
    call cross_product(b1,b2,b3)
    call vector_difference(p,v1,p1)
    call dotproduct(p1,b3,dp1b3)
    call dotproduct(b3,b3,db3b3)
    t3 = dp1b3/db3b3
    call scalar_product(t3,b3,t3b3)
    call vector_difference(p1,t3b3,q)
    call dotproduct(b1,b1,d11)
    call dotproduct(b1,b2,d12)
    call dotproduct(b2,b2,d22)
    call dotproduct(q,b1,q1)
    call dotproduct(q,b2,q2)
    t1 = (q1*d22 - q2*d12)/(d11*d22 - d12*d12)
    t2 = (d11*q2 - q1*d12)/(d11*d22 - d12*d12)
    if (t1 >= 0.0 .and. t2 >= 0.0 .and. (t1+t2) <= 1.0) then
       call pythag0(t3b3,d)
    else
       call minedge(v1,v2,p,d1)
       call minedge(v1,v3,p,d2)
       call minedge(v2,v3,p,d3)
       dx = amin1(d1,d2)
       d = amin1(dx,d3)
    end if

  end subroutine mintriangle

  !! ----------------------------------------------------------------

  !> Given a disk (center `c`, unit normal `n`, radius `r`), return
  !! the maximum distance from the origin to any point on the disk.
  !!
  !! @warning what does argument `t` stand for?
 
  subroutine maxdisk(c, n, r, t, dd)

    !!GCC$ ATTRIBUTES unused :: n

    !! Declare arguments
    real, dimension(3)   :: c, n    !> Disk center and normal
    real, intent(in)     :: r       !> Disk radius
    real, dimension(3,3) :: t       !> TBD
    real, intent(out)    :: dd      !> Maximum distance

    !! Local variables
    real, dimension(3) :: oi, q
    real               :: qq

    !! Get around compiler warning for unused dummy argument(s)
    real, dimension(3) :: n_2

    !! HACK !!
    if (.false.) then
       n_2(1) = n(1)
    end if

    call scalar_product(-1.0,c,oi)
    call rotate(q,t,oi)

    qq = sqrt(q(1)**2 + q(2)**2)
    dd = q(3)**2 + (qq + r)**2
    dd = sqrt(dd)

  end subroutine maxdisk

  !! ----------------------------------------------------------------

  !> @brief calculate the projection span of a disk onto an arbitary line
  !!
  !! @warning line passed through origin, returns two 1d quantities only.
  !!
  !! @warning what does argument `tr` stand for?

  subroutine prodisk(c, n, r, tr, v, s1, s2)

    !!GCC$ ATTRIBUTES unused :: n

    !! Declare arguments
    real, dimension(3)   :: c, n    !> Disk center and normal
    real, intent(in)     :: r       !> Disk radius
    real, dimension(3,3) :: tr      !> TBD
    real, dimension(3)   :: v       !> Line direction vector
    real, intent(out)    :: s1, s2  !> Two endpoints on line

    !! Local variables
    real, dimension(3) :: b, x1, x2
    real :: dot1, dot2, r1

    !! Get around compiler warning for unused dummy argument(s)
    real, dimension(3) :: n_2

    !! HACK !!
    if (.false.) then
       n_2(1) = n(1)
    end if

    call rotate(b,tr,v)

    r1 = sqrt(b(1)**2+b(2)**2)

    !> @warning compares two reals with ==
    if (r1 == 0.0) then
       x1(1) = r
       x1(2) = 0.0
       x1(3) = 0.0
    else
       x1(1) = r*b(1)/r1
       x1(2) = r*b(2)/r1
       x1(3) = 0.0
    end if
    x2(1) = -x1(1)
    x2(2) = -x1(2)
    x2(3) = 0.0
    call backtransform(c,tr,x1)
    call backtransform(c,tr,x2)
    call dotproduct(x1,v,dot1)
    call dotproduct(x2,v,dot2)
    s2 = amax1(dot1,dot2)
    s1 = amin1(dot1,dot2)

  end subroutine prodisk

  !! ----------------------------------------------------------------

  !> Given a disk (center `c`, unit normal `n`, radius `r`) and an
  !! arbitrary point `p`, compute the minimum distance between `p`
  !! and any point on the disk and return it in `d`.
  !!
  !! @warning what does argument `tr` stand for?
 
  subroutine mindisk(c, n, r, tr, p, d)

    !!GCC$ ATTRIBUTES unused :: tr

    !! Declare arguments
    real, dimension(3)   :: c, n     !> Disk center and unit normal
    real, intent(in)     :: r        !> Disk radius
    real, dimension(3,3) :: tr       !> TBD
    real, dimension(3)   :: p        !> Arbitrary point
    real, intent(out)    :: d        !> Minimum distance

    !! Local variables
    real, dimension(3) :: pmc, ta, nn
    real :: quid, t

    !! Get around compiler warning for unused dummy argument(s)
    real, dimension(3,3) :: tr_2

    !! HACK !!
    if (.false.) then
       tr_2(1,1) = tr(1,1)
    end if

    call vector_difference(p,c,pmc)
    call dotproduct(pmc,n,t)
    call scalar_product(t,n,ta)
    call vector_difference(pmc,ta,nn)
    call pythag0(nn,quid)

    if (quid <= r) then
       d = abs(t)
    else
       d = t**2 + (quid-r)**2
       d = sqrt(d)
    end if

  end subroutine mindisk

  !! ----------------------------------------------------------------
 
  !> Given a cylinder (center `c`, unit normal `n`, radius `r`,
  !! length `l`), compute the maximum distance of points on the
  !! surface of the cylinder to the origin.
  !!
  !! @warning what does argument `t` stand for?

  subroutine maxcylinder(c, n, r, l, t, dd)

    !!GCC$ ATTRIBUTES unused :: n

    !! Declare arguments
    real, dimension(3)   :: c, n     !> Cylinder center and unit normal
    real, intent(in)     :: r, l     !> Cylinder radius and length
    real, dimension(3,3) :: t        !> TBD
    real, intent(out)    :: dd       !> Maximum distance to origin

    !! Local variables
    real, dimension(3) :: oi, q
    real :: dd1, dd2, rr, zz

    !! Get around compiler warning for unused dummy argument(s)
    real, dimension(3) :: n_2

    !! HACK !!
    if (.false.) then
       n_2(1) = n(1)
    end if

    call scalar_product(-1.0,c,oi)
    call rotate(q,t,oi)

    rr = sqrt(q(1)**2 + q(2)**2)
    zz = q(3)
    dd1 = sqrt( (rr+r)**2 + (zz-l/2.0)**2 )
    dd2 = sqrt( (rr+r)**2 + (zz+l/2.0)**2 )
    dd = amax1(dd1,dd2)

  end subroutine maxcylinder

  !! ----------------------------------------------------------------

  !> Compute the projection span of the surface of a cylinder on an
  !! arbitrary line.
  !! 
  !! @warning what does argument `t` stand for?
  !!
  !! @warning arbitrary line goes through origin.

  subroutine procylinder(c, n, r, l, t, v, s1, s2)

    !! Declare arguments
    real, dimension(3)   :: c, n     !> Cylinder center and unit normal
    real, intent(in)     :: r, l     !> Cylinder radius and length
    real, dimension(3,3) :: t        !> TBD
    real, dimension(3)   :: v        !> Arbitrary line
    real, intent(out)    :: s1, s2   !> Projection span endpoints

    !! Local variables
    integer :: i
    real :: u1, u2, z1, z2
    real, dimension(3) :: c1, c2

    do i = 1,3
       c1(i) = c(i) + n(i)*l/2.0
       c2(i) = c(i) - n(i)*l/2.0
    end do

    call prodisk(c1,n,r,t,v,u1,u2)
    call prodisk(c2,n,r,t,v,z1,z2)
    s1 = amin1(u1,z1)
    s2 = amax1(u2,z2)

  end subroutine procylinder

  !! ----------------------------------------------------------------

  !> Compute the minimum distance between an arbitrary point `p` and
  !! points on the surface of a cylinder.
  !!
  !! @warning what does argument `t` stand for?

  subroutine mincylinder(c, n, r, l, t, p, d)

    !!GCC$ ATTRIBUTES unused :: n

    !! Declare arguments
    real, dimension(3)   :: c, n    !> Cylinder center and unit normal
    real, intent(in)     :: r, l    !> Cylinder radius and length
    real, dimension(3,3) :: t       !> TBD
    real, dimension(3)   :: p       !> Arbitrary point
    real, intent(out)    :: d       !> Minimum distance

    !! Local variables
    real, dimension(3) :: pmc, q
    real :: dp, ro, zz

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

    if (abs(q(3)) > l/2.0) then
       zz = abs(q(3)) - l/2.0
       d = sqrt(zz**2 + dp**2)
    else
       d = dp
    end if

  end subroutine mincylinder

  !! ----------------------------------------------------------------

  !> Compute the maximum distance between the origin and points on the
  !! surface of the torus.
  !!
  !! @verbatim
  !! TORUS[(cx,cy,cz),(nx,ny,nz),r1,r2]
  !! 
  !! Given a torus centered at c, with unit normal n, radii r1 and r2.
  !! This returns the maximum distance of points on the torus away 
  !! from the origin.
  !! 
  !!
  !!          x                         x
  !!      x       x                 x       x
  !!     x         x               x         x   ___
  !!     x         x               x         x    |
  !!      x       x                 x       x     |  r2
  !!          x                         x        ---
  !!
  !!          |---------2 * r1 ---------|
  !! @endverbatim
  !!
  !! @warning need a similar @e drawing for the other shapes.

  subroutine maxtorus(c, n, r1, r2, t, dd)

    !!GCC$ ATTRIBUTES unused :: n

    !! Declare arguments
    real, dimension(3) :: c, n      !> Torus center and unit normal
    real, intent(in)   :: r1, r2    !> Torus radii
    real, intent(out)  :: dd        !> Maximum distance to origin

    !! Local variables
    real, dimension(3) :: oi, q
    real, dimension(3,3) :: t
    real :: d1, d2, d3, d4, dod, sig, qq
    real :: x1, x2, y1, y2

    !! Get around compiler warning for unused dummy argument(s)
    real, dimension(3) :: n_2

    !! HACK !!
    if (.false.) then
       n_2(1) = n(1)
    end if

    call scalar_product(-1.0,c,oi)
    call rotate(q,t,oi)
    qq = sqrt(q(1)**2+q(2)**2)

    sig = q(3)/(qq-r1)
    dod = sqrt(r2**2/(1.0+sig**2))
    x1 = r1 + dod
    x2 = r1 - dod
    y1 = sig*dod
    y2 = -sig*dod
    d1 = (qq-x1)**2 + (q(3)-y1)**2
    d2 = (qq-x2)**2 + (q(3)-y2)**2

    sig = q(3)/(qq+r1)
    dod = sqrt(r2**2/(1.0+sig**2))
    x1 = -r1 + dod
    x2 = -r1 - dod
    y1 = sig*dod
    y2 = -sig*dod
    d3 = (qq-x1)**2 + (q(3)-y1)**2
    d4 = (qq-x2)**2 + (q(3)-y2)**2

    dd = amax1(d1,d2)
    dd = amax1(dd,d3)
    dd = amax1(dd,d4)

    dd = sqrt(dd)

  end subroutine maxtorus

  !! ----------------------------------------------------------------

  !> Compute the projection span of a torus onto an arbitrary line
  !!
  !! @warning What does argument `tr` stand for?
  !!
  !! @warning Line goes through the origin.

  subroutine protorus(c, n, r1, r2, tr, v, s1, s2)

    !!GCC$ ATTRIBUTES unused :: n

    !! Declare arguments
    real, dimension(3)   :: c, n       !> Torus center and unit normal
    real, intent(in)     :: r1, r2     !> Torus radii
    real, dimension(3,3) :: tr         !> TBD
    real, dimension(3)   :: v          !> Line direction vector
    real, intent(out)    :: s1, s2     !> Projection span endpoints

    real, dimension(3) :: b, x1, x2
    real :: cp, ct, dot1, dot2, sp, st

    !! Get around compiler warning for unused dummy argument(s)
    real, dimension(3) :: n_2

    !! HACK !!
    if (.false.) then
       n_2(1) = n(1)
    end if

    call rotate(b,tr,v)
    ct = b(3)
    st = sqrt(b(1)**2 + b(2)**2)

    !> @warning compares two reals with ==
    if (st == 0.0) then
       x1(1) = 0.0
       x1(2) = 0.0
       x1(3) = r2
    else
       cp = b(1)/st
       sp = b(2)/st
       x1(1) = (r1 + r2*st)*cp
       x1(2) = (r1 + r2*st)*sp
       x1(3) = r2*ct
    end if
    x2(1) = -x1(1)
    x2(2) = -x1(2)
    x2(3) = -x1(3)
    call backtransform(c,tr,x1)
    call backtransform(c,tr,x2)
    call dotproduct(x1,v,dot1)
    call dotproduct(x2,v,dot2)
    s2 = amax1(dot1,dot2)
    s1 = amin1(dot1,dot2)

  end subroutine protorus

  !! ----------------------------------------------------------------

  !> Compute the minimum distance between an arbitrary point and points
  !! on the surface of a torus.
  !! 
  !! @verbatim
  !! TORUS[(cx,cy,cz),(nx,ny,nz),r1,r2]
  !! 
  !! Given a torus centered at c, with unit normal n, radii r1 and r2.
  !! This returns the minimum distance of points on the torus away 
  !! from an arbitrary point p.
  !! 
  !! 
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
  !! 
  !! @warning what does argument `t` stands for?

  subroutine mintorus(c, n, r1, r2, t, p, d)

    !!GCC$ ATTRIBUTES unused :: n

    !! Declare arguments
    real, dimension(3)   :: c, n
    real, intent(in)     :: r1, r2
    real, dimension(3,3) :: t
    real, dimension(3)   :: p
    real, intent(out)    :: d

    !! Local variables
    real, dimension(3) :: pmc, q
    real :: d1, d2, d12, d22, qq

    !! Get around compiler warning for unused dummy argument(s)
    real, dimension(3) :: n_2

    !! HACK !!
    if (.false.) then
       n_2(1) = n(1)
    end if

    call vector_difference(p,c,pmc)
    call rotate(q,t,pmc)

    qq = sqrt(q(1)**2+q(2)**2)
    d12 = (qq-r1)**2 + q(3)**2
    d1 = sqrt(d12)
    d22 = (qq+r2)**2 + q(3)**2
    d2 = sqrt(d22)
    d = amin1(d1,d2)
    d = d - r2

  end subroutine mintorus

  !! ----------------------------------------------------------------

  !> @brief Compute maximum distance from origin to points on the cylinder.
  !! 
  !! Given a solid cylinder (center `c`, unit normal `n`, radius 
  !! `r`, length `l`), return the maximum distance of points on the
  !! cylinder away from the origin.
  !! 
  !! @warning what does argument `t` stand for?
  !! 
  !! @warning how does this routine differ from `maxcylinder`?
 
  subroutine maxsolcyl(c, n, r, l, t, dd)

    !!GCC$ ATTRIBUTES unused :: n

    !! Declare arguments
    real, dimension(3)   :: c, n     !> Cylinder center and unit normal
    real, intent(in)     :: r, l     !> Cylinder radius and length
    real, dimension(3,3) :: t        !> TBD
    real, intent(out)    :: dd       !> Maximum distance from origin

    !! Local variables
    real, dimension(3) :: oi, q
    real :: dd1, dd2, rr, zz

    !! Get around compiler warning for unused dummy argument(s)
    real, dimension(3) :: n_2

    !! HACK !!
    if (.false.) then
       n_2(1) = n(1)
    end if

    call scalar_product(-1.0,c,oi)
    call rotate(q,t,oi)

    rr = sqrt(q(1)**2 + q(2)**2)
    zz = q(3)
    dd1 = sqrt( (rr+r)**2 + (zz-l/2.0)**2 )
    dd2 = sqrt( (rr+r)**2 + (zz+l/2.0)**2 )
    dd = amax1(dd1,dd2)

  end subroutine maxsolcyl

  !! ----------------------------------------------------------------

  !> Compute the projection span of a solid cylinder onto an arbitary
  !! line
  !!
  !! @warning what does argument `t` stand for?
  !!
  !! @warning how does this routine differ from `procylinder`?
  !!
  !! @warning Assumes that line goes through the origin.

  subroutine prosolcyl(c, n, r, l, t, v, s1, s2)

    !! Declare arguments
    real, dimension(3)   :: c, n     !> Cylinder center and unit normal
    real, intent(in)     :: r, l     !> Cylinder radius and length
    real, dimension(3,3) :: t        !> TBD
    real, dimension(3)   :: v        !> Arbitrary line direction vector
    real, intent(out)    :: s1, s2   !> Projection span endpoints

    !! Local variables
    real, dimension(3) :: c1, c2
    real :: u1, u2, z1, z2
    integer :: i

    do i = 1,3
       c1(i) = c(i) + n(i)*l/2.0
       c2(i) = c(i) - n(i)*l/2.0
    end do

    call prodisk(c1,n,r,t,v,u1,u2)
    call prodisk(c2,n,r,t,v,z1,z2)
    s1 = amin1(u1,z1)
    s2 = amax1(u2,z2)

  end subroutine prosolcyl

  !! ----------------------------------------------------------------

  !> @brief Compute minimum distance between an arbitrary point and
  !! all points on a solid cylinder.
  !!
  !! Given a solid cylinder (center `c`, unit normal `n`, radius 
  !! `r`, length `l`) and an arbitrary point `p`, return the minimum
  !! distance between points on the cylinder and `p`.
  !! 
  !! @warning what does argument `t` stand for?
  !!
  !! @warning how does this routine differ from `mincylinder`?

  subroutine minsolcyl(c, n, r, l, t, p, d)

    !!GCC$ ATTRIBUTES unused :: n

    !! Declare arguments
    real, dimension(3)   :: c, n     !> Cylinder center and unit normal
    real, intent(in)     :: r, l     !> Cylinder radius and length
    real, dimension(3,3) :: t        !> TBD
    real, dimension(3)   :: p        !> Arbitrary point
    real, intent(out)    :: d        !> Minimum distance to `p`

    !! Local variables
    real, dimension(3) :: pmc, q
    real :: dp, ro, zz

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

    if (abs(q(3)) > l/2.0) then
       if (ro <= r) then
          d = abs(q(3)) - l/2.0
       else
          zz = abs(q(3)) - l/2.0
          d = sqrt(zz**2 + dp**2)
       end if
    else
       d = dp
    end if

  end subroutine minsolcyl

  !! ----------------------------------------------------------------

  !> Compute the maximum distance of points in the ellipsoid away from
  !! the origin.
  !! 
  !! @param dd -- return value.
  !!
  !! @warning how is the ellipsoid defined?
  !!
  !! @warning what does the argument `t` stand for?
  !!
  !! @warning explicitly inflate resut by 1.001.

  subroutine maxellipsoid(c, n1, n2, aa, bb, cc, t, dd)

    !! Declare arguments
    real, dimension(3)   :: c             !> Ellipsoid center
    real, dimension(3)   :: n1, n2        !> Ellipsoid two axes
    real, intent(in)     :: aa, bb, cc    !> Ellipsoid axial lengths
    real, dimension(3,3) :: t             !> TBD
    real, intent(out)    :: dd            !> Maximum distance to origin

    !! Local variables
    real, dimension(3) :: p
    integer :: minmax

    p(1) = 0.0
    p(2) = 0.0
    p(3) = 0.0

    minmax = 1
    call exellipsoid(c,n1,n2,aa,bb,cc,t,p,minmax,dd)

    !> Call to `exellipsoid` results in a value of `dd`
    !! that is always slightly low, so lets do this to
    !! compensate.  It only means we are using a slightly
    !! larger launch sphere.

    dd = dd*1.001

  end subroutine maxellipsoid

  !! ----------------------------------------------------------------

  !> Compute the projection span of an ellipsoid onto an arbitrary line.
  !!
  !! Upon entry, `v` is an arbitrary unit vector, upon which we are
  !! projecting the object
  !!
  !! @pre Upon entry, `c` is the center of the ellipsoid
  !!
  !! `t` is the rotation matrix required to rotate `n1` to x-axis,
  !! `n2` to y-axis, `n3` to z-axis

  subroutine proellipsoid(c, n1, n2, aa, bb, cc, t, v, s1, s2)

    !!GCC$ ATTRIBUTES unused :: n1, n2

    !! Declare arguments
    real, dimension(3)   :: c             !> Ellipsoid center
    real, dimension(3)   :: n1, n2        !> Ellipsoid two axes
    real, intent(in)     :: aa, bb, cc    !> Ellipsoid axial lengths
    real, dimension(3,3) :: t             !> Special rotation matrix
    real, dimension(3)   :: v             !> Direction vector of line
    real, intent(out)    :: s1, s2        !> Projection span endpoints

    !! Local variables
    real, dimension(3) :: vi, rp
    real ::    alam, zc
    integer :: i, j

    !! Get around compiler warning for unused dummy argument(s)
    real, dimension(3) :: n1_2, n2_2

    !! HACK !!
    if (.false.) then
       n1_2(1) = n1(1)
       n2_2(2) = n2(2)
    end if

    do i = 1,3
       vi(i) = 0.0
       do j = 1,3
          vi(i) = vi(i) + t(i,j)*v(j)
       end do
    end do

    alam = sqrt( (vi(1)*aa)**2 + (vi(2)*bb)**2 + (vi(3)*cc)**2 ) 
    alam = 1.0/alam

    rp(1) = alam*vi(1)*aa*aa
    rp(2) = alam*vi(2)*bb*bb
    rp(3) = alam*vi(3)*cc*cc

    s2 = rp(1)*vi(1) + rp(2)*vi(2) + rp(3)*vi(3) 
    s1 = -s2

    zc = c(1)*v(1) + c(2)*v(2) + c(3)*v(3)
    s1 = s1 + zc
    s2 = s2 + zc

  end subroutine proellipsoid

  !! ----------------------------------------------------------------

  !> Compute the minimum distance between the ellipsoid and a point `p`
  !!
  !! @warning return result in `d`.

  subroutine minellipsoid(c, n1, n2, aa, bb, cc, t, p, d)

    !! Declare arguments
    real, dimension(3)   :: c             !> Ellipsoid center
    real, dimension(3)   :: n1, n2        !> Ellipsoid two axes
    real, intent(in)     :: aa, bb, cc    !> Ellipsoid axial lengths
    real, dimension(3,3) :: t             !> Special rotation matrix
    real, dimension(3)   :: p             !> Arbitrary point
    real, intent(out)    :: d             !> Minimum distance

    !! Local variables
    integer :: minmax

    minmax = -1
    call exellipsoid(c,n1,n2,aa,bb,cc,t,p,minmax,d)

  end subroutine minellipsoid

  !! ----------------------------------------------------------------

  !> Given a line segment (`v1`, `v2`) and a point `p`, compute the
  !! minimum distance between the line segment and `p` and return it
  !! in `d`.

  subroutine minedge(v1, v2, p, d)

    !! Declare arguments
    real, dimension(3) :: v1, v2     !> Line segment vertices
    real, dimension(3) :: p          !> Arbitrary point
    real, intent(out)  :: d          !> Minimum distance

    !! Local variables
    real, dimension(3) :: p1, b, tb, n
    real :: b2, t, tb2

    call vector_difference(p,v1,p1)
    call vector_difference(v2,v1,b)
    call dotproduct(p1,b,tb2)
    call dotproduct(b,b,b2)
    t = tb2/b2

    if (t < 0.0) then
       call pythag(v1,p,d)
    else if (t > 1.0) then
       call pythag(v2,p,d)
    else
       call scalar_product(t,b,tb)
       call vector_difference(p1,tb,n)
       call pythag0(n,d)
    end if

  end subroutine minedge

  !! ================================================================

end module zeno_shape_distances

!!! ================================================================

!!! Local Variables:
!!! mode: f90
!!! time-stamp-line-limit: 30
!!! fill-column: 72
!!! End:
