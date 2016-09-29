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
!! This file provides the `zeno_integrations` module which defines
!! routines that perform ZENO integrations.
!!
!! @author Mark Mansfield (formerly of Stevens Institute) and students
!! @author Walid Keyrouz (walid.keyrouz@nist.gov)
!! @date   Fri Nov  1 13:53:09 2013 EDT
!!
!! @defgroup zeno_integrations Random-Walk Integrations
!! 
!! - Defines the following public entries that should be renamed: @c
!!   all_around, `wagaroo`, `blizzard`, `captain`, `calipers`.
!! 
!! - Defines the following private entries that may also be renamed
!!   later on to something that is more descriptive: `park`, 
!!   `greensphere`, `shine`, `reinit`.
!!
!! @todo Rename the routines in this module to have more meaningful
!! names.
!!
!! @warning The numbers `20` and `2000` seem to be too special!
!! 
! Time-stamp: <2015-01-05 16:43:58 walid>
! 
! ================================================================

module zeno_integrations

  !! ================================================================

  use numeric_kinds
  use numeric_constants

  use zeno_rng, only : get_rand, charge
  use zeno_stat_utils, only : mean20, tally, accume
  use zeno_shape_interior_test, only : inbody
  use zeno_shape_volume, only : primvol
  use zeno_shape_distances, only : distance, bridge
  use zeno_shape_surface_area, only : carea
  use zeno_shape_surface_point, only : getsurface
  use zeno_toss, only : toss_point
  use zeno_sphere, only : sphere
  use zeno_output, only : plus_or_minus

  implicit none
  private

  public all_around, wagaroo, blizzard, captain, calipers

contains

  !! ================================================================
  !!
  !! Public routines
  !!
  !! ================================================================

  !! @brief Do the interior integrations
  !! 
  !! Numerical integration of the radius of gyration and the volume
  !! 
  !! @warning Should rename routine to reflect functionality.
  !! @warning Argument list too long.

  subroutine all_around(maxelts, eltype, bv, nelts, &
       & xyzlow, xyzhih, rlaunch, rotations, m1, &
       & rg2int, delta_rg2int, volume, delta_volume, &
       & rg_done, mi, id, q, sq, dsq, round, savehits, tten)

    use zeno_debug, only : silent, print_flush
    use zeno_files_data

    !! Declare arguments
    integer, intent(in)          :: maxelts        !> Array capacity
    integer, dimension(maxelts)  :: eltype         !> Element types array
    real, dimension(maxelts,12)  :: bv             !> TBD
    integer, intent(in)          :: nelts          !> Array size
    real, dimension(3)           :: xyzlow         !> TBD
    real, dimension(3)           :: xyzhih         !> TBD
    real, intent(in)             :: rlaunch        !> Launch radius
    real, dimension(maxelts,3,3) :: rotations      !> Element rotations?
    integer, intent(in)          :: m1             !> Num of MC steps
    real, intent(out)            :: rg2int         !> TBD
    real, intent(out)            :: delta_rg2int   !> TBD
    real, intent(out)            :: volume         !> TBD
    real, intent(out)            :: delta_volume   !> TBD
    logical, intent(out)         :: rg_done        !> TBD
    integer, intent(out)         :: mi             !> TBD
    character(len=25)            :: id             !> Case ID
    real, dimension(0:81)        :: q              !> TBD
    real, dimension(0:81)        :: sq             !> TBD
    real, dimension(0:81)        :: dsq            !> TBD
    character(len=80)            :: round          !> TBD
    logical                      :: savehits       !> TBD
    real(dp_k), dimension(3,3)   :: tten           !> TBD

    !! Local variables
    real(dp_k), dimension(20) :: trials
    real(dp_k), dimension(20) :: rg2sum
    real(dp_k), dimension(20) :: rg2norm
    real(dp_k), dimension(0:81,20) :: sqsum
    real(dp_k) :: ttennorm
    real(dp_k) :: rij, qd, arg
    real, dimension(20) :: rad, sm1
    real, dimension(3)  :: rt1, rt2
    logical :: inside,early
    character(len=10)   :: mess
    character(len=2000) :: flush
    integer :: i, ii, i1, i2, j, jj, jj1, jj2, k, kk, loop
    integer :: mout, need
    real    :: an0, an1, ox, tn0, tn1, volbox

    if (savehits) open(unit=nih,file=fih,status='unknown')

    rg_done = .false.

    q(0) = 0.0
    do i = 1,81
       ox = float(i-41)/20.0
       q(i) = (10.0**ox)/rlaunch
    end do

    !> This is a stack of unprintable characters
    !! in order to flush the output buffer
    do i = 2,2000
       flush(i:i) = char(0)
    end do

    if (.not. silent) then
       write (*, '("INTERIOR CALCULATION"(4x)(a25)(i12))') id, m1
       write (*, '(72("="))')
    end if

    mout = 0

    do k = 1,20
       trials(k) = 0.0d0
       sm1(k) = 0.0d0
       rg2sum(k) = 0.0d0
       rg2norm(k) = 0.0d0
       do ii = 0,81
          sqsum(ii,k) = 0.0d0
       end do
    end do

    do jj1 = 1,3
       do jj2 = 1,3
          tten(jj1,jj2) = 0.0d0
       end do
    end do
    ttennorm = 0.0d0

    Do_MC_Steps: do i = 1,m1
       loop = mod(i,20) + 1
       sm1(loop) = sm1(loop) + 1.0

       inside = .false.
       do while (.not. inside)
          do ii = 1,3
             rt1(ii) = xyzlow(ii) + get_rand()*(xyzhih(ii)-xyzlow(ii))
          end do

          !> rt1 is the trial point

          call inbody(rt1,maxelts,eltype,bv,nelts,rotations,inside,early,mess)

          if (early) then
             write(nzno, '((a), (a10), (a))') &
                  & "  Body element of type: ", mess, &
                  & "  is inconsistent with INTERIOR calculation."
             if (savehits) close(nih)
             return
          end if
          
          trials(loop) = trials(loop) + 1.0d0
       end do

       sm1(loop) = sm1(loop) + 1.0

       inside = .false.
       do while (.not. inside)
          do ii = 1,3
             rt2(ii) = xyzlow(ii) + get_rand()*(xyzhih(ii)-xyzlow(ii))
          end do

          !> rt2 is the second trial point

          call inbody(rt2,maxelts,eltype,bv,nelts,rotations,inside,early,mess)

          trials(loop) = trials(loop) + 1.0d0
       end do

       !> At this point, we have two points, rt1 and rt2, and we know
       !! that both of them lie inside the body.
       
       if (savehits) then
          write(nih, '(3(g20.8))') rt1
          write(nih, '(3(g20.8))') rt2
       end if

       rij = 0.0d0
       do jj = 1,3
          rij = rij + dble(rt1(jj)-rt2(jj))**2
          rg2sum(loop) = rg2sum(loop) + dble(rt1(jj)-rt2(jj))**2
       end do
       rg2norm(loop) = rg2norm(loop) + 2.0d0
       do jj1 = 1,3
          do jj2 = 1,3
             tten(jj1,jj2) = tten(jj1,jj2) + &
                  &          (rt1(jj1)-rt2(jj1))*(rt1(jj2)-rt2(jj2))
          end do
       end do

       ttennorm = ttennorm + 2.0d0
       rij = dsqrt(rij)
       do kk = 0,81
          qd = dble(q(kk))
          arg = qd*rij
          if (kk == 0) then
             sqsum(kk,loop) = sqsum(kk,loop) + 1.0d0
          else
             sqsum(kk,loop) = sqsum(kk,loop) + dsin(arg)/arg
          end if
       end do

       need = nint(79.0*float(i)/float(m1))

       do while (mout < need)
          flush(1:1) = round(mout+1:mout+1)

          if (.not. silent .and. print_flush) then
             write(*,"((a2000))",advance="no") flush
          end if

          mout = mout + 1
       end do

    end do Do_MC_Steps

    if (.not. silent) write(*,"( )")

    volbox = 1.0
    do i = 1,3
       volbox = volbox *(xyzhih(i)-xyzlow(i))
    end do

    do i = 1,20
       rad(i) = volbox*sngl(sm1(i)/trials(i))
    end do
    call mean20(rad,volume,delta_volume)

    do i = 1,20
       rad(i) = sngl(rg2sum(i)/rg2norm(i))
    end do
    call mean20(rad,rg2int,delta_rg2int)

    do i1 = 1,3     
       do i2 = 1,3   
          tten(i1,i2) = tten(i1,i2)/ttennorm
       end do
    end do

    do i = 1,20
       rad(i) = sngl(sqsum(0,i))
    end do

    call mean20(rad,an0,an1)
    sq(0) = 1.0
    dsq(0) = an1/an0

    do j = 1,81
       do i = 1,20
          rad(i) = sngl(sqsum(j,i))
       end do
       call mean20(rad,tn0,tn1)
       sq(j) = tn0/an0
       dsq(j) = tn1/an0
    end do

    rg_done = .true.
    mi = m1
    
    if (savehits) close(nih)

  end subroutine all_around

  !! ----------------------------------------------------------------

  !> @brief Do the covered-interior integration
  !! 
  !! Numerical integration of the radius of gyration and the volume
  !! 
  !! @warning Definitely change the name.
  !! @warning Argument list too long.

  subroutine wagaroo(maxelts, eltype, bv, nelts, rlaunch, &
       &             rotations, m1, rg2int, delta_rg2int, &
       &             volume, delta_volume, rg_done, mi, id, &
       &             q, sq, dsq, round, savehits, vaar, tten)

    use zeno_debug, only : silent
    use zeno_files_data

    !! Declare arguments
    integer, intent(in)          :: maxelts        !> Array capacity
    integer, dimension(maxelts)  :: eltype         !> Element types array
    real, dimension(maxelts,12)  :: bv             !> TBD
    integer, intent(in)          :: nelts          !> Array size
    real, intent(in)             :: rlaunch        !> Launch radius
    real, dimension(maxelts,3,3) :: rotations      !> Element rotations?
    integer, intent(in)          :: m1             !> Num of MC steps
    real, intent(out)            :: rg2int         !> TBD
    real, intent(out)            :: delta_rg2int   !> TBD
    real, intent(out)            :: volume         !> TBD
    real, intent(out)            :: delta_volume   !> TBD
    logical, intent(out)         :: rg_done        !> TBD
    integer, intent(out)         :: mi             !> TBD
    character(len=25)            :: id             !> Case ID
    real, dimension(0:81)        :: q              !> TBD
    real, dimension(0:81)        :: sq             !> TBD
    real, dimension(0:81)        :: dsq            !> TBD
    character(len=80)            :: round          !> TBD
    logical, intent(in)          :: savehits       !> TBD
    real, dimension(maxelts)     :: vaar           !> TBD
    real(dp_k), dimension(3,3)   :: tten           !> Global rotation matrix?

    !! Local variables
    real(dp_k), dimension(20) :: trials
    real(dp_k), dimension(20) :: rg2sum, rg2norm
    real(dp_k), dimension(0:81,20) :: sqsum
    real(dp_k) :: ttennorm
    real(dp_k) :: rij,qd,arg
    real, dimension(12)  :: pass_vec
    real, dimension(3,3) ::tumble
    real, dimension(20)  :: rad
    real, dimension(3)   :: rt1, rt2
    character(len=10)    :: mess
    character(len=2000)  :: flush
    logical :: early
    integer :: i, ii, ii1, ii2, j, jj, jj1, jj2, k, kk, kk1, kk2
    integer :: loop, mout, myelt, need
    real    :: an0, an1, ox, result, tn0, tn1

    volume = 0.0
    delta_volume = 0.0

    Do_Nelts: do i = 1,nelts

       myelt = eltype(i)
       do j = 1,12
          pass_vec(j) = bv(i,j)
       end do
       do j = 1,3
          do k = 1,3
             tumble(j,k) = rotations(i,j,k)
          end do
       end do
       call primvol(myelt,pass_vec,tumble,result,early,mess)
       vaar(i) = result
       volume = volume + result

       if (early) then
          write(nzno, '("  Body element of type: ",(a10))') mess
          write(nzno, '("  is inconsistent with C1-INTERIOR calculation.")')

          return
       end if

    end do Do_Nelts
          
    if (savehits) open(unit=nih,file=fih,status='unknown')

    rg_done = .false.

    q(0) = 0.0
    do i = 1,81
       ox = float(i-41)/20.0
       q(i) = (10.0**ox)/rlaunch
    end do

    !> This is a stack of unprintable characters
    !! in order to flush the output buffer
    do i = 2,2000
       flush(i:i) = char(0)
    end do

    if (.not. silent) then
       write(*, "(('C1 INTERIOR CALCULATION'), (1x), (a25), (i12))") id, m1
       write(*, "(79('='))")
    end if

    mout = 0

    do k = 1,20
       trials(k) = 0.0d0
       rg2sum(k) = 0.0d0
       rg2norm(k) = 0.0d0
       do ii = 0,81
          sqsum(ii,k) = 0.0d0
       end do
    end do

    do kk1 = 1,3    
       do kk2 = 1,3
          tten(kk1,kk2) = 0.0d0
       end do
    end do
    ttennorm = 0.0d0

    Do_MC_Steps: do i = 1,m1

       loop = mod(i,20) + 1

       !> rt1 is the first point
       call toss_point(rt1,maxelts,eltype,bv,nelts,rotations,vaar,volume)
       trials(loop) = trials(loop) + 1.0d0

       !> rt2 is the second point
       call toss_point(rt2,maxelts,eltype,bv,nelts,rotations,vaar,volume)

       trials(loop) = trials(loop) + 1.0d0

       !> At this point, we have two points, rt1 and rt2, and we know
       !! that both of them lie inside the body.

       if (savehits) then
          write(nih, '(3(g20.8))') rt1
          write(nih, '(3(g20.8))') rt2
       end if

       rij = 0.0d0
       do jj = 1,3
          rij = rij + dble(rt1(jj)-rt2(jj))**2
          rg2sum(loop) = rg2sum(loop) + dble(rt1(jj)-rt2(jj))**2
       end do

       rg2norm(loop) = rg2norm(loop) + 2.0d0
       do jj1 = 1,3
          do jj2 = 1,3
             tten(jj1,jj2) = tten(jj1,jj2) + &
                  &          (rt1(jj1)-rt2(jj1))*(rt1(jj2)-rt2(jj2))
          end do
       end do

       ttennorm = ttennorm + 2.0d0
       rij = dsqrt(rij)

       do kk = 0,81
          qd = dble(q(kk))
          arg = qd*rij
          if (kk == 0) then
             sqsum(kk,loop) = sqsum(kk,loop) + 1.0d0
          else
             sqsum(kk,loop) = sqsum(kk,loop) + dsin(arg)/arg
          end if
       end do

       need = nint(79.0*float(i)/float(m1))
       do while (mout < need)
          flush(1:1) = round(mout+1:mout+1)
          if (.not. silent) write(*, '((a2000))', advance="no") flush
          mout = mout + 1
       end do

    end do Do_MC_Steps

    if (.not. silent) write(*,'( )')

    do i = 1,20
       rad(i) = sngl(rg2sum(i)/rg2norm(i))
    end do
    call mean20(rad,rg2int,delta_rg2int)

    do ii1 = 1,3
       do ii2 = 1,3
          tten(ii1,ii2) = tten(ii1,ii2)/ttennorm
       end do
    end do

    do i = 1,20
       rad(i) = sngl(sqsum(0,i))
    end do

    call mean20(rad,an0,an1)
    sq(0) = 1.0
    dsq(0) = an1/an0

    do j = 1,81
       do i = 1,20
          rad(i) = sngl(sqsum(j,i))
       end do
       call mean20(rad,tn0,tn1)
       sq(j) = tn0/an0
       dsq(j) = tn1/an0
    end do

    rg_done = .true.
    mi = m1

    if (savehits) close(nih)

  end subroutine wagaroo

  !! ----------------------------------------------------------------

  !> @brief Does ZENO integrations.
  !!
  !! Generate random walk trajectories to do path-integration
  !!
  !! @warning How does this differ from the other integration routines?
  !! 
  !! @warning Should be renamed to reflect functionality.
  !! 
  !! @warning Uses anonymous numeric constants.  `nneb` argument set to
  !! `3` before routine is called.
  !! 
  !! @warning Too many parameters

  subroutine blizzard(maxelts, eltype, bv, nelts, m1, tol, rlaunch, &
       &              rotations, cap, alpha_bongo, tol_given, zeno_done, &
       &              delta_cap, delta_bongo, mz, id, tae, uae, vae, wae, &
       &              round, bubble, bubble_rad, nebtab, nneb, ninn, &
       &              rlist,  strikes, savehits)

    use zeno_debug, only : silent, print_flush
    use zeno_files_data

    !! Declare arguments
    integer, intent(in)            :: maxelts     !> Array capacity
    integer, dimension(maxelts)    :: eltype      !> ELements array
    real, dimension(maxelts,12)    :: bv          !> TBD
    integer, intent(in)            :: nelts       !> Array size
    integer, intent(in)            :: m1          !> Num of MC steps
    real, intent(in out)           :: tol         !> TBD
    real, intent(in)               :: rlaunch     !> Launch radius
    real, dimension(maxelts,3,3)   :: rotations   !> Element rotation matrices
    real                           :: cap         !> TBD
    real, dimension(3,3)           :: alpha_bongo !> TBD
    logical                        :: tol_given   !> TBD
    logical                        :: zeno_done   !> TBD
    real                           :: delta_cap   !> TBD
    real, dimension(3,3)           :: delta_bongo !> TBD
    integer                        :: mz          !> TBD
    character(len=25)              :: id          !> Case ID
    real(dp_k)                     :: tae         !> TBD
    real(dp_k), dimension(3)       :: uae         !> TBD
    real(dp_k), dimension(3,3)     :: vae         !> TBD
    real(dp_k), dimension(3,3)     :: wae         !> TBD
    character(len=80)              :: round       !> TBD
    real, dimension(3)             :: bubble      !> TBD
    real                           :: bubble_rad  !> TBD
    integer, dimension(maxelts)    :: nebtab      !> TBD
    integer                        :: nneb        !> TBD
    integer                        :: ninn        !> TBD
    real, dimension(maxelts)       :: rlist       !> TBD
    real(dp_k), dimension(maxelts) :: strikes     !> TBD
    logical                        :: savehits    !> TBD
        
    !! Local variables
    real(dp_k), dimension(3,3,20) :: vp, ve
    real(dp_k), dimension(20)     :: sum
    real(dp_k), dimension(3)      :: taer
    real(dp_k) :: grand
    real, dimension(3,3) :: copy, dopy
    real, dimension(3,3) :: aa, daa
    real, dimension(3)   :: rt
    integer, dimension(3,20) :: khitp, khite
    integer, dimension(3)    :: kk
    integer :: hitelt
    logical :: hit
    character(len=2000) :: flush
    character(len=3)    :: pom
    integer :: i, ii, j, jj, jax, loop, mout, mtdo, need
    real    :: r, r2

    if (savehits) open(unit=nzh,file=fzh,status='unknown')

    zeno_done = .false.

    if (.not. tol_given) then
       tol = rlaunch/1.0e6
       tol_given = .true.
    end if

    r = rlaunch
    r2 = r**2
    mtdo = m1

    !> This is a stack of unprintable characters
    !! in order to flush the output buffer
    do i = 2,2000
       flush(i:i) = char(0)
    end do

    do i = 1,nelts
       strikes(i) = 0.0d0
    end do

    call reinit(khitp,khite,vp,ve,sum)

    if (.not. silent) then
       write(*, '(("ZENO CALCULATION"), (8x), (a25), (i12))') id, m1
       write(*, '(79("="))')
    end if

    mout = 0

    Do_MC_Steps: do jax = 1,mtdo

       loop = mod(jax,20) + 1
       sum(loop) = sum(loop) + 1.0d0

       call sphere(rt,r)
       call charge(rt,r,kk)
       call park(maxelts,eltype,bv,rotations,nelts,rt,r,r2,hit,tol, &
            &    bubble,bubble_rad,nebtab,nneb,ninn,rlist,hitelt)

       if (hit) then
          call accume(rt,kk,khitp,khite,vp,ve,loop)
          strikes(hitelt) = strikes(hitelt) + 1.0d0
          if (savehits) then
             call plus_or_minus(kk,pom)
             write(nzh, '((a3), 3(g20.8))') pom, rt
          end if
       end if

       need = nint(79.0*float(jax)/float(mtdo))
       do while (mout < need)
          flush(1:1) = round(mout+1:mout+1)

          if (.not. silent .and. print_flush) then
             write(*, '((a2000))', advance='no') flush
          end if

          mout = mout + 1
       end do

    end do Do_MC_Steps

    if (.not. silent) write(*, '( )')

    call tally(khitp,khite,vp,ve,sum,aa,daa,cap,delta_cap,r, &
         &     taer,uae,vae,wae,rlaunch)

    tae = taer(1)
    do i = 1,3
       do j = 1,3
          wae(i,j) = wae(i,j) * (12.0d0 * M_PI) * (dble(r)**2)
          vae(i,j) = vae(i,j) * (12.0d0 * M_PI) * (dble(r)**2)
       end do
    end do

    !> Because of sampling error, the polarizability tensor is not
    !! exactly symmetric.  Now we symmetrize it.

    do ii = 1,3
       do jj = 1,3
          copy(ii,jj) = aa(ii,jj)
          dopy(ii,jj) = daa(ii,jj)
       end do
    end do

    aa(1,2) = 0.5*(aa(1,2)+aa(2,1))
    aa(1,3) = 0.5*(aa(1,3)+aa(3,1))
    aa(2,3) = 0.5*(aa(2,3)+aa(3,2))
    aa(2,1) = aa(1,2)
    aa(3,1) = aa(1,3)
    aa(3,2) = aa(2,3)
    daa(1,2) = 0.5*(daa(1,2)+daa(2,1))
    daa(1,3) = 0.5*(daa(1,3)+daa(3,1))
    daa(2,3) = 0.5*(daa(2,3)+daa(3,2))
    daa(2,1) = daa(1,2)
    daa(3,1) = daa(1,3)
    daa(3,2) = daa(2,3)

    do ii = 1,3
       do jj = 1,3
          alpha_bongo(ii,jj) = aa(ii,jj)
          delta_bongo(ii,jj) = daa(ii,jj)
       end do
    end do

    zeno_done = .true.
    mz = mtdo

    open(unit=nstk,file=fstk,status='unknown')
    grand = 0.0d0
    do i = 1,nelts
       grand = grand + strikes(i)
       write(nstk, '((i10), (f20.0))') i, strikes(i)
    end do

    write(nstk,"('Grand total:  ',(f20.0))") grand

    close(nstk)

    if (savehits) close(nzh)

  end subroutine blizzard

  !! ----------------------------------------------------------------

  !> @brief Perform surface integrations.
  !! 
  !! Perform surface integrations by computing Kirkwood estimates of the
  !! hydrodynamic radius and its uncertainty.  The name is `captain` as
  !! an allusion to Startrek's Captain Kirk.
  !!
  !! @warning argument list is too long.
  !!
  !! @warning should be renamed to something more meaninful.

  subroutine captain(maxelts, eltype, bv, nelts, m1do, rotations, &
       &             kirk_done, saar, kirk, delta_kirk, &
       &             surface, delta_surface, rg2surf, delta_rg2surf, &
       &             ms, id, round, savehits)

    use zeno_debug, only : silent
    use zeno_codes_data
    use zeno_files_data

    !> Declare arguments
    integer, intent(in)          :: maxelts        !> Array capacity
    integer, dimension(maxelts)  :: eltype         !> Elements array
    real, dimension(maxelts,12)  :: bv             !> Elts body values
    integer, intent(in)          :: nelts          !> Array size
    integer, intent(in)          :: m1do           !> Num of MC steps
    real, dimension(maxelts,3,3) :: rotations      !> Elt rotation matrices
    logical, intent(out)         :: kirk_done      !> Surf. calculation flag
    real, dimension(maxelts)     :: saar           !> Surf. areas
    real, intent(out)            :: kirk           !> Kirkwood hydrodyn. rad.
    real, intent(out)            :: delta_kirk     !> Uncertainty in `kirk`,
    real, intent(out)            :: surface        !> Surface area
    real, intent(out)            :: delta_surface  !> Surface area uncertainty
    real, intent(out)            :: rg2surf        !> Square radius gyration
    real, intent(out)            :: delta_rg2surf  !> Sq. rad. gur. uncertain.
    integer, intent(out)         :: ms             !> TBD
    character(len=25)            :: id             !> Case ID
    character(len=80)            :: round          !> TBD
    logical, intent(in)          :: savehits       !> CLI option flag?

    !! Local variables

    real(dp_k), dimension(20) :: sum1, sum2, trials, successes
    real(dp_k), dimension(20) :: rg2sum, rg2norm
    real,       dimension(3)  :: v1, v2
    real,       dimension(20) :: rad
    real    :: delta_rho, rho, rr, rr2, total
    integer :: i, j, k, loop, mout, ncub, need, npil

    character(len=2000) :: flush
            
    if (savehits) open(unit=nsh,file=fsh,status='unknown')

    kirk_done = .false.

    npil = 0
    ncub = 0
    do i = 1,nelts
       if (eltype(i) == pillar_code) npil = npil + 1
       if (eltype(i) == cube_code) ncub = ncub + 1
    end do

    if (npil > 1) then
       write(nzno,"('SURFACE integration may be unreliable whenever there')")
       write(nzno,"('are abutting pillars.')")
    end if

    if (ncub > 1) then
       write(nzno,"('SURFACE integration may be unreliable whenever there')")
       write(nzno,"('are abutting cubes.')")
    end if

    call carea(maxelts,eltype,bv,nelts,saar,total)
    do i = 1,20
       sum1(i) = 0.0d0
       sum2(i) = 0.0d0
       trials(i) = 0.0d0
       successes(i) = 0.0d0
       rg2sum(i) = 0.0d0
       rg2norm(i) = 0.0d0
    end do

    !> This is a stack of unprintable characters
    !> in order to flush the output buffer
    do i = 2,2000
       flush(i:i) = char(0)
    end do

    if (.not. silent) then
       write(*, "('SURFACE CALCULATION',(5x),(a25),(i12))") id, m1do
       write(*, "(79('='))")
    end if

    mout = 0

    Do_MC_Steps: do i = 1,m1do

       loop = mod(i,20) + 1
       call getsurface(maxelts,eltype,bv,nelts,saar,total,v1, &
            &          trials,rotations,loop)
       call getsurface(maxelts,eltype,bv,nelts,saar,total,v2, &
            &          trials,rotations,loop)
       if (savehits) then
          write(nsh,'(3(g20.8))') v1
          write(nsh,'(3(g20.8))') v2
       end if

       successes(loop) = successes(loop) + 2.0d0
       rr2 = 0.0
       do j = 1,3
          rr2 = rr2 + (v1(j)-v2(j))**2
       end do
       rr = sqrt(rr2)
       sum2(loop) = sum2(loop) + 1.0d0/dble(rr)
       sum1(loop) = sum1(loop) + 1.0d0
       rg2sum(loop) = rg2sum(loop) + dble(rr2)
       rg2norm(loop) = rg2norm(loop) + 2.0d0

       need = nint(79.0*float(i)/float(m1do))
       do while (mout < need)
          flush(1:1) = round(mout+1:mout+1)
          if (.not. silent) write(*,'(a2000)',advance='no') flush
          mout = mout + 1
       end do

    end do Do_MC_Steps

    if (.not. silent ) write(*, '(" ")')

    do k = 1,20
       rad(k) = total*sngl(successes(k)/trials(k))
    end do
    call mean20(rad,surface,delta_surface)

    do i = 1,20
       rad(i) = sngl(rg2sum(i)/rg2norm(i))
    end do
    call mean20(rad,rg2surf,delta_rg2surf)

    do k = 1,20
       rad(k) = sngl(sum2(k)/sum1(k))
    end do
    call mean20(rad,rho,delta_rho)
    kirk = 1.0/rho
    delta_kirk = delta_rho/(rho**2)

    kirk_done = .true.
    ms = m1do

    if (savehits) close(nsh)

  end subroutine captain

  !! ----------------------------------------------------------------

  !> @brief Do the project-onto-line integration.
  !!
  !! @warning Should be renamed!
  !!
  !! @warning Argument list is too long!
  !!
  !! @warning uses anonymous numeric constants.

  subroutine calipers(maxelts, eltype, bv, nelts, m1, rotations, &
       &              span_done, span, delta_span, &
       &              shadow_done, shadow, delta_shadow, &
       &              rlaunch, mp, id, round, savehits)

    use zeno_debug, only : silent
    use zeno_codes_data
    use zeno_files_data

    !! Declare arguments
    integer, intent(in)          :: maxelts      !> Array capacity
    integer, dimension(maxelts)  :: eltype       !> Elements array
    real, dimension(maxelts,12)  :: bv           !> Elt. body values array
    integer, intent(in)          :: nelts        !> Array size
    integer, intent(out)         :: m1           !> Num of MC steps
    real, dimension(maxelts,3,3) :: rotations    !> Elt rotations array
    logical, intent(out)         :: span_done    !> Proj. on line calc flag
    real, intent(out)            :: span         !> Giddings length
    real, intent(out)            :: delta_span   !> Giddigns length uncertainty
    logical, intent(out)         :: shadow_done  !> Proj. on surf. calc flag 
    real, intent(out)            :: shadow       !> TBD
    real, intent(out)            :: delta_shadow !> TBD
    real, intent(in)             :: rlaunch      !> Launch radius
    integer, intent(out)         :: mp           !> TBD
    character(len=25)            :: id           !> Case ID
    character(len=80)            :: round        !> TBD
    logical, intent(in)          :: savehits     !> CLI option flag?

    !! Local variables
    real(dp_k), dimension(20) :: sum, acd, shoo
    real,       dimension(20) :: x 
    real,       dimension(3)  :: v
    character(len=2000)       :: flush

    logical :: do_shadow
    integer :: i, jax, kin, loop, mout, mtdo, need
    real    :: d, rlsh

    if (savehits) open(unit=nph,file=fph,status='unknown')

    span_done = .false.
    shadow_done = .false.
    rlsh = M_PI_sp * rlaunch * rlaunch

    !> Attempt to calculate planar projection only if the
    !> body is composed completely of spheres

    do_shadow = .true.
    do i = 1,nelts
       if (eltype(i) /= sphere_code) then
          do_shadow = .false.
          exit
       end if
    end do

    do i = 2,2000
       flush(i:i) = char(0)
    end do

    do i = 1,20
       sum(i) = 0.0d0
       acd(i) = 0.0d0
       shoo(i) = 0.0d0
    end do

    if (.not. silent) then
       write(*,"('PROJECTION CALCULATION',(2x),(a25),(i12))") id,m1
       write(*,"(79('='))")
    end if

    mout = 0
    mtdo = m1

    Do_MC_Steps: do jax = 1,mtdo

       loop = mod(jax,20) + 1
       sum(loop) = sum(loop) + 1.0d0
       call sphere(v,1.0)
       call bridge(maxelts,eltype,bv,nelts,rotations,v,d)

       if (savehits) write(nph,'(g20.8)') d

       acd(loop) = acd(loop) + dble(d)
       if (do_shadow) then
          call shine(maxelts,bv,nelts,v,rlaunch,kin)
          shoo(loop) = shoo(loop) + dble(kin)
       end if
       need = nint(79.0 * float(Jax)/float(mtdo))
       do while (mout < need)
          flush(1:1) = round(mout+1:mout+1)

          if (.not. silent) write(*,'((a2000))',advance='no') flush

          mout = mout + 1
       end do

    end do Do_MC_Steps

    if (.not. silent) write(*,'( )')

    do i = 1,20
       x(i) = sngl(acd(i)/sum(i))
    end do
    call mean20(x,span,delta_span)
    span_done = .true.

    if (do_shadow) then
       do i = 1,20
          x(i) = sngl(shoo(i)/sum(i))
       end do
       call mean20(x,shadow,delta_shadow)
       shadow_done = .true.
       shadow = shadow * rlsh
       delta_shadow = delta_shadow * rlsh
    end if

    mp = mtdo

    if (savehits) close(nph)

  end subroutine calipers

  !! ================================================================
  !!
  !! Private routines
  !!
  !! ================================================================

  !> @brief subroutine `park`
  !!
  !! Upon entry, the random walker sits at a point `p` on the launch
  !! sphere.  `r` is the radius of the launch sphere, and `r2` is its
  !! square.
  !!
  !! This subroutine lets the walker drift down onto the object or else
  !! out to infinity.  If it is lost to infinity, `hit` is returned as
  !! `.false.`, and if it drifts onto the body, `hit` is returned as
  !! `.true.` and `p` is returned with the point at which the walker
  !! hits.
  !!
  !! This is the algorithm obeyed by the subroutine:
  !! 
  !! -#
  !!   If the point `p` lies inside the launch sphere, proceed to step
  !!   2.  If it lies outside the launch sphere, move it onto the launch
  !!   sphere using the charge-outside-a-sphere Green's function, which
  !!   may also move the point off to infinity.  If the point gets moved
  !!   off to infinity, set `hit = .false.` and return.
  !!
  !! -#
  !!   The walker is now on or inside the launch sphere, but outside the
  !!   body.  A call to distance returns the distance to the body, 
  !!   `ds`.  If `ds` is less than `tol`, the program assumes that the
  !!   walker has adsorbed and so `hit` is set equal to `.true.` and
  !!   we return.  Otherwise, we jump to the surface of the sphere that
  !!   is centered on the current point and has radius `ds`.
  !!
  !! -#
  !!   Loop back to step 1.
  !!
  !! @warning Argument list is too long!

  subroutine park(maxelts, eltype, bv, rotations, nelts, p, &
       &          r, r2, hit, tol, bubble, bubble_rad, &
       &          nebtab, nneb, ninn, rlist, hitelt)

    !! Declare arguments
    integer, intent(in)          :: maxelts
    integer, dimension(maxelts)  :: eltype
    real, dimension(maxelts,12)  :: bv
    real, dimension(maxelts,3,3) :: rotations
    integer, intent(in)          :: nelts
    real, dimension(3)           :: p
    real, intent(in)             :: r, r2
    logical, intent(out)         :: hit
    real, intent(in)             :: tol
    real, dimension(3)           :: bubble
    real, intent(out)            :: bubble_rad
    integer, dimension(maxelts)  :: nebtab
    integer                      :: nneb
    integer                      :: ninn
    real, dimension(maxelts)     :: rlist
    integer, intent(out)         :: hitelt

    !! Local variables
    real, dimension(3) :: d
    real               :: ds, r0
    logical :: gone
    integer :: i, nearto
    integer :: pass

    !> On the first time, a new bubble
    !> is always needed
    pass = 0
    bubble_rad = -1.0

    Passes_Lbl: do

       pass = pass + 1

       !> STEP 1:

       r0 = p(1)**2 + p(2)**2 + p(3)**2

       if (r0 > r2) then

          !> No need to call greensphere, on first pass
          !! it could only be  due to round-off that we are here.

          if (pass == 1) then
             gone = .false.
          else
             r0 = sqrt(r0)
             call greensphere(p,r,r0,gone)

             !> A new bubble is needed everytime
             !! we reposition the launch sphere
             bubble_rad = -1.0
          end if

          if (gone) then
             hit = .false.
             return
          end if

       end if                   ! if (r0 > r2) then

       !> STEP 2:

       call distance(maxelts, eltype, bv, nelts, rotations, p, ds, &
            &        bubble, bubble_rad, nebtab, nneb, ninn, rlist, nearto)

       if (ds < tol) then
          hit = .true.        
          hitelt = nearto
          return
       end if

       call sphere(d,ds)
       do i = 1,3
          p(i) = p(i) + d(i)     !  Tentative new position
       end do

    end do Passes_Lbl

  end subroutine park

  !! ----------------------------------------------------------------

  !> @brief Subroutine `greensphere`
  !!
  !! @warning to be renamed
  !!
  !! Upon entry, `p` contains the coordinate of the point, `r`
  !! contains the radius of the launch sphere, and `r0` contains the
  !! distance of `p` from the origin.  It is known that `p` lies
  !! outside the launch sphere.
  !!
  !! This determines the probability that the walker escapes to infinity
  !! without ever returning to the launch sphere, and lets the walker
  !! escape with that probablity.  Otherwise, it returns the walker to a
  !! new point on the launch sphere.
  !!
  !! Upon exit, `gone = .true.` if the walker has escaped to
  !! infinity.  Otherwise, `gone = .false.`, and `p` contains the
  !! coordinates of the new point on the launch sphere.

  subroutine greensphere(p, r, r0, gone)

    !! Declare arguments
    real, dimension(3)   :: p
    real, intent(in)     :: r, r0
    logical, intent(out) :: gone

    !! Local variables
    integer :: i
    real    :: alpha, cosp, cost, cp, cx, phir, pip
    real    :: s, sinp, sint, sp, sx
    real    :: t1, t2, t3, thetr, xn, yn, zn

    alpha = r/r0
    gone = get_rand() > alpha
    if (gone) return

    do i = 1,3
       p(i) = p(i)/r0
    end do

    !> Pull off the numbers needed to transform the p-vector to the
    !! plus z-axis:

    cost = p(3)

    !> Minor correction for round-off errors:
    if (cost > 1.0) cost = 1.0
    if (cost < -1.0) cost = -1.0

    thetr = acos(cost)
    phir = atan2(p(2),p(1))

    !> With modern fortran implementations, atan2 has a range of 2*pi,
    !! so further correction of phir is not necessary

    !> For the time being, assume a new coordinate system for which
    !! the radius vector is on the z-axis (thetr and phir will permit
    !! us to transform back later)

    s = get_rand()

    if (s == 0.0) then
       cx = -1.0
    else
       t1 = (1.0 + alpha**2)/(2.0*alpha)
       t2 = (1.0 - alpha**2)**2
       t3 = (1.0 - alpha + 2.0*alpha*s)**2
       cx = t1 - t2/(2.0*alpha*t3)
    end if

    !> Minor correction for roundoff errors:
    if (cx > +1.0) cx = +1.0
    if (cx < -1.0) cx = -1.0

    sx = sqrt(1.0 - cx**2)
    pip = 2.0 * M_PI_sp * get_rand()
    sp = sin(pip)
    cp = cos(pip)

    !> This is the new position of the diffusor in 
    !! the 2nd coordinate system

    p(1) = r*sx*cp
    p(2) = r*sx*sp
    p(3) = r*cx

    !> Now transform to the first coordinate system, first
    !! by a rotation about y-axis through thetr

    sint = sin(thetr)
    zn = p(3)*cost - p(1)*sint
    xn = p(3)*sint + p(1)*cost
    p(1) = xn
    p(3) = zn

    !> and then by a rotation about z-axis through phir

    cosp = cos(phir)
    sinp = sin(phir)

    xn = p(1)*cosp - p(2)*sinp
    yn = p(2)*cosp + p(1)*sinp
    p(1) = xn
    p(2) = yn

  end subroutine greensphere

  !! ----------------------------------------------------------------

  !> @brief Subroutine `shine` TBD
  !! 
  !! @warning uses anonymous numeric constants.

  subroutine shine(maxelts, bv, nelts, v, rlaunch, kin)

    !! Declare arguments
    integer, intent(in)         ::  maxelts
    real, dimension(maxelts,12) :: bv
    integer, intent(in)         :: nelts
    real, dimension(3)          :: v
    real, intent(in)            :: rlaunch
    integer, intent(out)        :: kin

    !! Local variables
    real, dimension(3) :: vn, vp, vm, vs
    real, dimension(3) :: c, pc
    integer :: i, j
    real    :: alpha, beta, rad, rpick, rt2, tpick, x, y, z

    kin = 0

    !> construct an arbitrary vn, normal to v
    call sphere(vp,1.0)
    beta = 0.0
    do j = 1,3      
       beta = beta + v(j)*vp(j)
    end do
    do j = 1,3
       vn(j) = vp(j) - beta*v(j)
    end do
    beta = vn(1)*vn(1) + vn(2)*vn(2) + vn(3)*vn(3)
    beta = sqrt(beta)
    vn(1) = vn(1)/beta 
    vn(2) = vn(2)/beta 
    vn(3) = vn(3)/beta 

    !> construct vm, normal to both v and vn
    vm(1) = v(2)*vn(3) - v(3)*vn(2)
    vm(2) = v(3)*vn(1) - v(1)*vn(3)
    vm(3) = v(1)*vn(2) - v(2)*vn(1)

    !> Pick a point at random inside the projection of the launch
    !! radius

    z = get_rand()
    rpick = rlaunch*sqrt(z)
    tpick = 2.0 * M_PI_sp * get_rand()
    x = rpick*cos(tpick)
    y = rpick*sin(tpick)

    do i = 1,3
       vs(i) = x*vn(i) + y*vm(i)
    end do

    Do_Nelts : do i = 1,nelts

       do j = 1,3
          c(j) = bv(i,j)
       end do
       rad = bv(i,4)
       alpha = 0.0
       do j = 1,3
          alpha = alpha + v(j)*c(j)
       end do
       do j = 1,3
          pc(j) = c(j) - alpha*v(j)
       end do

       rt2 = 0.0
       do j = 1,3
          rt2 = rt2 + (pc(j)-vs(j))**2
       end do
       rt2 = sqrt(rt2)
       if (rt2 < rad) then
          kin = 1
          return
       end if

    end do Do_Nelts

  end subroutine shine

  !! ----------------------------------------------------------------

  !> @brief initialize statistical registers
  !! 
  !! @warning should be renamed!

  subroutine reinit(khitp, khite, vp, ve, sum)

    !! Declare arguments
    integer, dimension(3,20)      :: khitp    !> TBD
    integer, dimension(3,20)      :: khite    !> TDB
    real(dp_k), dimension(3,3,20) :: vp       !> TBD
    real(dp_k), dimension(3,3,20) :: ve       !> TBD
    real(dp_k), dimension(20)     :: sum      !> TBD

    !! Local variables
    integer :: i, j, k

    do k = 1,20
       do i = 1,3
          do j = 1,3
             vp(i,j,k) = 0.0d0
             ve(i,j,k) = 0.0d0
          end do
          khite(i,k) = 0
          khitp(i,k) = 0
       end do
       sum(k) = 0.0d0
    end do

  end subroutine reinit

  !! ================================================================

end module zeno_integrations

!!! ================================================================

!!! Local Variables:
!!! mode: f90
!!! time-stamp-line-limit: 30
!!! fill-column: 72
!!! End:
