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
!! This file provides the `zeno_ellipsoid_utils` module which defines
!! routines related to the `exellipsoid` shape routine in the
!! file @c "zeno-shape-dist".
!!
!! @author Mark Mansfield (formerly of Stevens Institute) and students
!! @author Walid Keyrouz (walid.keyrouz@nist.gov)
!! @date   Wed Oct 23 18:11:05 2013 EDT
!! 
!! @defgroup zeno_ellipsoid_utils Ellipsoid Computation Utils
!! 
!! Groups together routines used by the `exellipsoid` routine in the
!! file "zeno-shape-dist".  As of the time of creation of this module,
!! these routines are only called by the `exellipsoid` routine.  An
!! alternative is to simply have these routines as `private` entries in
!! the "zeno-shape-dist" module.
!! 
!! @warning this module's routines will most likely undergo dramatic
!! change later on.
!! 
!! @warning this module's routines should be redefined as part of a
!! class hierarchy.
!! 
! Time-stamp: <2015-01-05 16:43:43 walid>
! 
! ================================================================

module zeno_ellipsoid_utils

  !! ================================================================

  use numeric_kinds

  use zeno_vectors, only : vector_difference, rotate
  use zeno_sort_utils, only : sort

  implicit none
  !! private sort, above, between, eval, eval2, converge, deval
  private

  public exellipsoid

contains

  !! ================================================================
  !!
  !! Public routines
  !!
  !! ================================================================
 
  !! @brief determines the extreme points on the surface of an ellipsoid.
  !! 
  !! The determination of extremal points on the surface of an
  !! ellipsoid `p` is difficult when any component of `p` is close
  !! to zero, because the problem becomes degenerate.
  !! 
  !! We deal with the problem here by making any very small
  !! coordinates not so small.  The justification for this adjustment
  !! is different in each of two different cases:
  !! @enumerate
  !! @item All three coordinates near zero, and
  !! @item Only one or two coordinates near zero.
  !! @endenumerate
  !! 
  !! When all three coordinates are near zero, we change them to a
  !! small positive value epsilon.  This increases the maximum
  !! distance by an amount of order epsilon.  We are obviously at a
  !! point inside the ellipsoid, and so we are attempting to set the
  !! radius of the launch sphere.  Therefore, the only effect of this
  !! manipulation is to make the launch sphere slightly larger.
  !! 
  !! If we are not on the interior of the ellipsoid trying to set the
  !! launch sphere radius, then we must be on the exterior and trying
  !! to compute the minimum distance to the surface from an exterior
  !! point.  Therefore, the second case above only occurs when we are
  !! on the outside looking in.  Then, if we change only one or two
  !! zero coordinates by a small amount `epsilon`, the effect on the
  !! distance measurement is of order @f$\epsilon^2$@f (it's like we are
  !! rotating a lever arm) and therefore negligible.
  !! 
  !! @warning routine may be renamed.
  !! @warning too many arguments.

  subroutine exellipsoid(c, n1, n2, aa, bb, cc, t, po, minmax, d)

    !!GCC$ ATTRIBUTES unused :: n1, n2

    !! Declare arguments
    real, dimension(3)   :: c
    real, dimension(3)   :: n1, n2
    real, intent(in)     :: aa, bb, cc
    real, dimension(3,3) :: t
    real, dimension(3)   :: po
    integer, intent(in)  :: minmax
    real, intent(out)    :: d

    !! Local variables
    real, dimension(3) :: pomc
    real, dimension(3) :: a, a2s
    real, dimension(3) :: p, x, x1, x2
    real, dimension(6,3) :: xsol
    integer :: i, j, ndif, ns, nsol
    real ::    alam, alam1, alam2, dd, ss

    !! Get around compiler warning for unused dummy argument(s)
    real, dimension(3) :: n1_2, n2_2

    !! HACK !!
    if (.false.) then
       n1_2(1) = n1(1)
       n2_2(2) = n2(2)
    end if

    call vector_difference(po,c,pomc)
    call rotate(p,t,pomc)

    a(1) = aa
    a(2) = bb
    a(3) = cc
        
    !> The determination of extremal points on the surface of an
    !> ellipsoid is difficult when any component of p is close to
    !> zero, because the problem becomes degenerate.
    
    !> We deal with the problem here by making any very small
    !> coordinates not so small.  The justification for this
    !> adjustment is different in each of two different cases: 1.  All
    !> three coordinates near zero, and 2. Only one or two coordinates
    !> near zero.
    
    !> When all three coordinates are near zero, we change them to a
    !> small positive value epsilon.  This increases the maximum
    !> distance by an amount of order epsilon.  We are obviously at a
    !> point inside the ellipsoid, and so we are attempting to set the
    !> radius of the launch sphere.  Therefore, the only effect of
    !> this manipulation is to make the launch sphere slightly larger.
    
    !> If we are not on the interior of the ellipsoid trying to set
    !> the launch sphere radius, then we must be on the exterior and
    !> trying to compute the minimum distance to the surface from an
    !> exterior point. Therefore, the second case above only occurs
    !> when we are on the outside looking in.  Then, if we change only
    !> one or two zero coordinates by a small amount epsilon, the
    !> effect on the distance measurement is of order epsilon**2 (it's
    !> like we are rotating a lever arm) and therefore negligible.
    
    do i = 1,3
       if (abs(p(i))/a(i) < 1.0e-6) p(i) = 1.0e-4 * a(i)
    end do

    nsol = 0
    !> a2s contains squared and sorted a's:
    !> ndif = 3:  a2s(1) < a2s(2) < a2s(3)
    !> ndif = 2:  a2s(1) < a2s(2) and a2s(3) is irrelevant
    !> ndif = 1:  sphere, a2s(1) is only relevant entry
    call sort(a,a2s,ndif)

    if (ndif == 3) then
            
       call above(a,p,alam,x,-a2s(1),+1.0)
       nsol = nsol + 1
       xsol(nsol,1) = x(1)
       xsol(nsol,2) = x(2)
       xsol(nsol,3) = x(3)

       call between(a,p,alam1,alam2,x1,x2,-a2s(2),-a2s(1),ns)
       if (ns == 2) then
          xsol(nsol+1,1) = x1(1)
          xsol(nsol+1,2) = x1(2)
          xsol(nsol+1,3) = x1(3)
          xsol(nsol+2,1) = x2(1)
          xsol(nsol+2,2) = x2(2)
          xsol(nsol+2,3) = x2(3)
          nsol = nsol + 2
       end if

       call between(a,p,alam1,alam2,x1,x2,-a2s(3),-a2s(2),ns)
       if (ns == 2) then
          xsol(nsol+1,1) = x1(1)
          xsol(nsol+1,2) = x1(2)
          xsol(nsol+1,3) = x1(3)
          xsol(nsol+2,1) = x2(1)
          xsol(nsol+2,2) = x2(2)
          xsol(nsol+2,3) = x2(3)
          nsol = nsol + 2
       end if

       call above(a,p,alam,x,-a2s(3),-1.0)
       nsol = nsol + 1
       xsol(nsol,1) = x(1)
       xsol(nsol,2) = x(2)
       xsol(nsol,3) = x(3)

    else if (ndif == 2) then
                
       call above(a,p,alam,x,-a2s(1),+1.0)
       nsol = nsol + 1
       xsol(nsol,1) = x(1)
       xsol(nsol,2) = x(2)
       xsol(nsol,3) = x(3)

       call between(a,p,alam1,alam2,x1,x2,-a2s(2),-a2s(1),ns)
       if (ns == 2) then
          xsol(nsol+1,1) = x1(1)
          xsol(nsol+1,2) = x1(2)
          xsol(nsol+1,3) = x1(3)
          xsol(nsol+2,1) = x2(1)
          xsol(nsol+2,2) = x2(2)
          xsol(nsol+2,3) = x2(3)
          nsol = nsol + 2
       end if

       call above(a,p,alam,x,-a2s(2),-1.0)
       nsol = nsol + 1
       xsol(nsol,1) = x(1)
       xsol(nsol,2) = x(2)
       xsol(nsol,3) = x(3)

    else

       call above(a,p,alam,x,-a2s(1),+1.0)
       nsol = nsol + 1
       xsol(nsol,1) = x(1)
       xsol(nsol,2) = x(2)
       xsol(nsol,3) = x(3)

       call above(a,p,alam,x,-a2s(1),-1.0)
       nsol = nsol + 1
       xsol(nsol,1) = x(1)
       xsol(nsol,2) = x(2)
       xsol(nsol,3) = x(3)

    end if

    do i =  1, nsol
       dd = 0.0
       ss = 0.0
       do j = 1,3
          dd = dd + (p(j) - xsol(i,j))**2
          ss = ss + (xsol(i,j)/a(j))**2
       end do
       dd = sqrt(dd)
       if (i == 1) then
          d = dd
       else
          if (minmax == -1) then
             d = amin1(dd,d)
          else
             d = amax1(dd,d)
          end if
       end if
    end do

  end subroutine exellipsoid

  !! ================================================================
  !!
  !! Private routines
  !!
  !! ================================================================

  !! @brief TBD

  subroutine above(a, p, alam, x, sing, ccc)

    !! Declare arguments
    real, dimension(3) :: a
    real, dimension(3) :: p
    real, dimension(3) :: x
    real               :: alam
    real               :: sing
    real               :: ccc

    !! Local variables
    real :: add, t, z, zhih, zlow, ztry, zz

    add = 1.0

    do
       z = sing + add*ccc
       call eval(z,a,p,x,t)
       if (t < 1.0) then
          zhih = z
          exit
       end if
       add = add * 1.5
    end do

    ztry = zhih

    do
       ztry = (ztry+sing)/2.0
       call eval(ztry,a,p,x,t)
       if (t > 1.0) then
          zlow = ztry
          exit
       end if
    end do

    if (ccc < 0.0) then
       zz = zlow
       zlow = zhih
       zhih = zz
    end if

    call converge(zlow,zhih,a,p,x,alam)

  end subroutine above

  !! ----------------------------------------------------------------

  !> @brief TBD
  !! 
  !! @warning Numeric constants used directly in code!

  subroutine between(a, p, alam1, alam2, x1, x2, sing1, sing2, ns)

    !! Declare arguments
    real, dimension(3) :: a, p
    real               :: alam1, alam2
    real, dimension(3) :: x1, x2
    real               :: sing1, sing2
    integer            :: ns

    !! Local variables
    real, dimension(3) :: x
    real    :: t, z, z1, zc, zhih, zlow, zmid, ztest
    integer :: k1

    z = (sing1+sing2)/2.0
    call eval2(z,a,p,x,t)
        
    if (t > 0.0) then
       zlow = z
       k1 = 1
    else
       zhih = z        
       k1 = 2
    end if

    if (k1 == 1) then

       do
          ztest = (zlow+sing2)/2.0
          call eval2(ztest,a,p,x,t)
          if (t < 0.0) then
             zhih = ztest
             exit
          end if
          zlow = ztest
       end do

    else

       do
          ztest = (zhih + sing1)/2.0
          call eval2(ztest,a,p,x,t)
          if (t > 0.0) then
             zlow = ztest
             exit
          end if
          zhih = ztest
       end do

    end if                      ! if (k1 == 1) then

    do while ((zhih-zlow)*2.0/(abs(zhih)+abs(zlow)) > 1.0e-5)

       zmid = (zhih + zlow)/2.0
       call eval2(zmid,a,p,x,t)
       if (t < 0.0) then
          zhih = zmid
       else
          zlow = zmid
       end if

    end do

    zc = 0.5*(zhih + zlow)
    call eval(zc,a,p,x,t)
    if (t > 1.0) then
       ns = 0
       return
    end if

    ns = 2
    z1 = zc

    do
       call eval(z1,a,p,x,t)
       if (t > 1.0) then
          exit
       end if
       z1 = (z1 + sing1)/2.0
    end do

    call converge(z1,zc,a,p,x1,alam1)

    z1 = zc

    do
       call eval(z1,a,p,x,t)
       if (t > 1.0) then    
          exit
       end if
       z1 = (z1 + sing2)/2.0
    end do

    call converge(zc,z1,a,p,x2,alam2)

  end subroutine between

  !! ----------------------------------------------------------------

  !> @brief TBD
  !>
  !> @warning may be renamed to something more meaningful

  subroutine eval(z, a, p, x, t)

    !! Declare arguments
    real               :: z
    real, dimension(3) :: a, p, x
    real               :: t

    !! Local variables
    integer :: i
    real    :: q

    do i = 1,3
       q = 1.0 + z/(a(i)**2)
       q = 1.0/q
       x(i) = q*p(i)
    end do
        
    t = 0.0
    do i = 1,3
       t = t + (x(i)/a(i))**2
    end do

  end subroutine eval

  !! ----------------------------------------------------------------

  !> @brief `double` version of `eval`
  !>
  !> @warning may be renamed to something more meaningful

  subroutine deval(z, a, p, x, t)

    !! Declare arguments

    real(dp_k)               :: z
    real(dp_k), dimension(3) :: a, p, x
    real(dp_k)               :: t

    !! Local variables
    real(dp_k) :: q
    integer  :: i

    do i = 1,3
       q = 1.0d0 + z/(a(i)**2)
       q = 1.0d0/q
       x(i) = q*p(i)
    end do
        
    t = 0.0d0
    do i = 1,3
       t = t + (x(i)/a(i))**2
    end do

  end subroutine deval

  !! ----------------------------------------------------------------

  !> Subroutine `eval2`
  !>
  !> @brief TBD
  !>
  !> @warning what are the differences between `eval`, `eval2`, and `deval`?
  !>
  !> @warning may be renamed to something more meaningful

  subroutine eval2(z, a, p, x, t)

    !! Declare arguments
    real :: z
    real, dimension(3) :: a, p, x
    real :: t

    !! Local variables
    integer :: i
    real    :: den, q

    do i = 1,3
       q = 1.0 + z/(a(i)**2) 
       q = 1.0/q
       x(i) = q*p(i)
    end do

    t = 0.0 
    do i = 1,3
       den = (a(i)**2 + z)**2
       t = t + x(i)*p(i)/den
    end do

  end subroutine eval2

  !! ----------------------------------------------------------------

  !> @brief TBD
  !>
  !> @warning uses anonymous numeric constants in code!

  subroutine converge(zlow, zhih, a, p, x, alam)

    !! Declare arguments
    real :: zlow, zhih
    real, dimension(3) :: a, p, x
    real :: alam

    !! Local variables
    real(dp_k), dimension(3) :: da, dp, dx
    real(dp_k) :: dhih, dlow, dm, dt, dthih, dtlow
    integer  :: i
    real     :: t

    dhih = dble(zhih)
    dlow = dble(zlow)

    do i = 1,3
       da(i) = dble(a(i))
       dp(i) = dble(p(i))
    end do

    call deval(dhih,da,dp,dx,dthih)
    call deval(dlow,da,dp,dx,dtlow)

    do while ((dhih-dlow)*2.0d0/(dabs(dhih)+dabs(dlow)) > 1.0d-12)
       dm = 0.5d0*(dhih + dlow)
       call deval(dm,da,dp,dx,dt)
       if ( (dthih-1.0d0)*(dt-1.0d0) > 0.0d0) then
          dhih = dm
          dthih = dt
       else
          dlow = dm
          dtlow = dt
       end if
    end do

    zhih = sngl(dhih)
    zlow = sngl(dlow)

    alam = (zhih + zlow)/2.0        
    call eval(alam,a,p,x,t)

  end subroutine converge

  !! ================================================================

end module zeno_ellipsoid_utils

!!! ================================================================

!!! Local Variables:
!!! mode: f90
!!! time-stamp-line-limit: 30
!!! fill-column: 72
!!! End:
