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
!! This file defines the `zeno_output` module which provides output
!! functionality for ZENO.
!! 
!! @author Mark Mansfield (formerly of Stevens Institute) and students
!! @author Walid Keyrouz (walid.keyrouz@nist.gov)
!! @date   Mon Sep 23 14:16:36 2013 EDT
!! 
!! @defgroup zeno_output Output Functionality
!! 
!! Groups together ZENO's output routines and defines two entry points:
!! `begin_output`, `report`, `alldone`, and `assemble`. and 
!! `plus_or_minus`.  All other routines are internal.
!!
!! Private routines are: show_errors, poplar, writeno, dotout,
!! writeyes, writeme, longjohn, pack20, pack49, shiftleft, 
!! shiftright, derive, bundle.
!!
!! @warning Added hacks to get around compiler warnings about unused
!! dummy arguments in subroutines `report` and `polar`
!!
!! @warning Uses too many literal constants to denote string sizes and
!! others.
!! 
! Time-stamp: <2015-01-05 16:44:20 walid>
!
! ================================================================

module zeno_output

  !! ================================================================

  use numeric_kinds
  use numeric_constants

  use zeno_math_utils
  use zeno_kernels, only : compvisc, lighthalf, loadconstants, &
       & rudnick_gaspari
  use zeno_units, only : conventional, derive, diffuse, friction, &
       & parameters, svedberg
  use zeno_tensors, only : ttdiag

  implicit none
  private
  public begin_output, report, alldone, gettime, assemble, &
       & plus_or_minus

contains

  !! ================================================================
  !!
  !! Public routines
  !!
  !! ================================================================

  !> Sends initial information to report file (`nzno` I/O unit)
  !! 
  !! @param id    (in) ZENO case identifier
  !! @param start (in) Starting date and time (Unix style)
  !! @param nelts (in) Number of body elements

  subroutine begin_output(id, nelts, start)

    use zeno_files_data

    !! Declare arguments

    character(len=25) :: id     !> ZENO case identifier
    character(len=28) :: start  !> Starting date and time (Unix style)
    integer           :: nelts  !> Number of body elements

    write(nzno,"('START (Version 4.0): ',(a28))") start
    write(*,   "('START (Version 4.0): ',(a28))") start

    write(nzno,"('Body name:  ',(a25))") id
    write(nzno,"('Number of body elements: ',(i7))") nelts
    write(nzno,"(50('='))")

  end subroutine begin_output

  !! ----------------------------------------------------------------

  !> @brief Generates program output to report file.
  !!
  !! Generates program output to report file.  Assumes that preliminary
  !! data was already written by `begin_output`.  This subroutine
  !! writes everything else.
  !!
  !! @pre assumes that `begin_output` has already been called.
  !!
  !! @warning too many parameters!
  !! @warning may need to add more boolean parameters to control output
  !!
  !! @warning compiler generates two unused dummy arguments: id, nelts
  !!
  !! The following parameters are booleans and store information
  !! about the computation and control reporting of results.
  !!
  !!  - `bt` -- Temperature was set in body file
  !!  - `bm` -- Mass was set in body file
  !!  - `bw` -- Solvent=water command was found in body file
  !!  - `bc` -- Solvent viscosity was set in body file
  !!  - `bbf` -- Buoyancy factor was set in body file

  subroutine report(id, actions, m1, nelts, tol, rlaunch, &
       &            cap, delta_cap, alpha_bongo, delta_bongo, tten, &
       &            volume, delta_volume, surface, delta_surface, &
       &            rg2int, delta_rg2int, rg2surf, delta_rg2surf, &
       &            kirk, delta_kirk, span, delta_span, &
       &            shadow, delta_shadow, shadow_done, &
       &            zeno_done, kirk_done, rg_done, span_done, &
       &            launch_done, tol_given, unitcode, mz, mi, ms, mp, &
       &            bt, bm, bw, bc, bbf, temp, tunit, mass, &
       &            munit, visc, vunit, buoy, hscale, &
       &            tae, uae, vae, wae, q2pade, eigens, xx, &
       &            q, sq, dsq, xyzlow, xyzhih)

    !!GCC$ ATTRIBUTES unused :: id, nelts

    use zeno_files_data

    !! Declare arguments
    character(len=25)        :: id             !> case identifier
    character(len=4)         :: actions        !> action specifications
    integer, dimension(4)    :: m1             !> TBD
    integer                  :: nelts          !> number of bodies
    real                     :: tol            !> TBD
    real                     :: rlaunch        !> TBD
    real                     :: cap            !> TBD
    real                     :: delta_cap      !> TBD
    real, dimension(3,3)     :: alpha_bongo    !> TBD
    real, dimension(3,3)     :: delta_bongo    !> TBD
    real(dp_k), dimension(3,3) :: tten         !> TBD
    real                     :: volume         !> TBD
    real                     :: delta_volume   !> TBD
    real                     :: surface        !> TBD
    real                     :: delta_surface  !> TBD
    real                     :: rg2int         !> TBD
    real                     :: delta_rg2int   !> TBD
    real                     :: rg2surf        !> TBD
    real                     :: delta_rg2surf  !> TBD
    real                     :: kirk           !> TBD
    real                     :: delta_kirk     !> TBD
    real                     :: span           !> TBD
    real                     :: delta_span     !> TBD
    real                     :: shadow         !> TBD
    real                     :: delta_shadow   !> TBD
    logical                  :: shadow_done    !> TBD
    logical                  :: zeno_done      !> TBD
    logical                  :: kirk_done      !> TBD
    logical                  :: rg_done        !> TBD
    logical                  :: span_done      !> TBD
    logical                  :: launch_done    !> TBD
    logical                  :: tol_given      !> TBD
    character(len=2)         :: unitcode       !> TBD
    integer                  :: mz             !> TBD
    integer                  :: mi             !> TBD
    integer                  :: ms             !> TBD
    integer                  :: mp             !> TBD
    logical                  :: bt             !> Temp set in body file
    logical                  :: bm             !> Mass set in body file
    logical                  :: bw             !> Solvent=water in body file
    logical                  :: bc             !> Solvent viscosity in body file
    logical                  :: bbf            !> Buoancy factor in body file
    real(dp_k), dimension(2) :: temp           !> TBD
    character(len=6)         :: tunit          !> TBD
    real(dp_k), dimension(2) :: mass           !> TBD
    character(len=6)         :: munit          !> TBD
    real(dp_k), dimension(2) :: visc           !> TBD
    character(len=6)         :: vunit          !> TBD
    real(dp_k), dimension(2) :: buoy           !> TBD
    real                     :: hscale         !> TBD
    real(dp_k)               :: tae            !> TBD
    real(dp_k), dimension(3) :: uae            !> TBD
    real(dp_k), dimension(3,3) :: vae          !> TBD
    real(dp_k), dimension(3,3) :: wae          !> TBD
    real                     :: q2pade         !> TBD
    real(dp_k), dimension(3) :: eigens         !> TBD
    real(dp_k), dimension(2) :: xx             !> TBD
    real, dimension(0:81)    :: q              !> TBD
    real, dimension(0:81)    :: sq             !> TBD
    real, dimension(0:81)    :: dsq            !> TBD
    real, dimension(3)       :: xyzlow         !> TBD
    real, dimension(3)       :: xyzhih         !> TBD

    !! Local variables

    !> The following local boolean variables complement the boolean
    !! parameters listed above and store information about the
    !! computation and control reporting of results.
    !!
    !! - `bs` -- Surface integration succeeded
    !! - `bi` -- Interior integration succeeded
    !! - `bl` -- Launch radius was calculated
    !! - `bk` -- Skin thickness was determined
    !! - `bz` -- Zeno integration succeeded
    !! - `bp` -- Projection-onto-line integration succeeded
    !! - `bpp` -- Projection-onto-plane integration succeeded

    logical :: bs, bi, bl, bk, bz, bp, bpp

    !! Additional local variables

    real(dp_k), dimension(2) :: capd
    character(len=6)         :: lunit,qunit
    real, dimension(3,3)     :: alpha_paper, delta_paper
    real                     :: pi
    character(len=1)         :: b1,b2
    logical                  :: do_something
    real(dp_k), dimension(3) :: box1, box2
    real(dp_k)               :: rl, eps
    real(dp_k), dimension(2) :: a11, a12, a13
    real(dp_k), dimension(2) :: a21, a22, a23
    real(dp_k), dimension(2) :: a31, a32, a33
    real(dp_k), dimension(2) :: rk, surf, rg2surfd
    real(dp_k), dimension(2) :: v, rg2intd
    real(dp_k), dimension(2) :: lhalf
    real(dp_k), dimension(2) :: qa, qb, sqa, sqb
    real(dp_k), dimension(2) :: rp, rap
    real(dp_k)               :: h1, h2, h3
    real(dp_k), dimension(0:81,2) :: stfac
    integer                  :: i, i1, j, ka, kb, nzip
    logical                  :: l_half_flg

    !! Get around compiler warnings for unused arguments
    character(len=25) :: id_2
    integer :: nelts_2

    !! HACK!!!
    if (.false.) then
       id_2 = id
       nelts_2 = nelts
    end if

    pi = M_PI_sp
    h1 = dble(hscale)
    h2 = h1**2
    h3 = h1**3

    lunit(1:2) = unitcode
    lunit(3:6) = '    '      !  Length unit
    qunit = lunit
    qunit(3:6) = '^-1 '

    call show_errors()          ! Report detected errors in output file

    !! Send a job summary to the output file

    write(nzno,"('JOB SUMMARY:')")
    write(nzno,"('       Actions      Checked if    Monte Carlo')")
    write(nzno,"('       requested    completed     size')")
    write(nzno,"(50('-'))")

    !> Output format
    !! @verbatim
    !!        Actions      Checked if    Monte Carlo
    !!        requested    completed     size
    !! ==================================================
    !!            z            *          1000000
    !! 000000000111111111122222222223333333333444
    !! 123456789012345678901234567890123456789012
    !! @endverbatim

    do i = 1,4
       do_something = .false.
       if (actions(i:i) == 'z') then
          b1 = 'z'
          b2 = ' '
          if (zeno_done) b2 = '*'
          nzip = m1(i)
          do_something = .true.
       else if (actions(i:i) == 'Z') then
          b1 = 'Z'
          b2 = ' '
          if (zeno_done) b2 = '*'
          nzip = m1(i)
          do_something = .true.
       else if (actions(i:i) == 'i') then
          b1 = 'i'
          b2 = ' '
          if (rg_done) b2 = '*'
          nzip = m1(i)
          do_something = .true.
       else if (actions(i:i) == 'I') then
          b1 = 'I'
          b2 = ' '
          if (rg_done) b2 = '*'
          nzip = m1(i)
          do_something = .true.
       else if (actions(i:i) == 'c') then
          b1 = 'c'
          b2 = ' '
          if (rg_done) b2 = '*'
          nzip = m1(i)
          do_something = .true.
       else if (actions(i:i) == 'C') then
          b1 = 'C'
          b2 = ' '
          if (rg_done) b2 = '*'
          nzip = m1(i)
          do_something = .true.
       else if (actions(i:i) == 's') then
          b1 = 's'
          b2 = ' '
          if (kirk_done) b2 = '*'
          nzip = m1(i)
          do_something = .true.
       else if (actions(i:i) == 'S') then
          b1 = 'S'
          b2 = ' '
          if (kirk_done) b2 = '*'
          nzip = m1(i)
          do_something = .true.
       else if (actions(i:i) == 'p') then
          b1 = 'p'
          b2 = ' '
          if (span_done) b2 = '*'
          nzip = m1(i)
          do_something = .true.
       else if (actions(i:i) == 'P') then
          b1 = 'P'
          b2 = ' '
          if (span_done) b2 = '*'
          nzip = m1(i)
          do_something = .true.
       end if

       if (do_something) then
          write(nzno,"((11x),(a1),(12x),(a1),(7x),(i10))") b1, b2, nzip
       end if

    end do

    write(nzno,"(50('='))")

    bl = launch_done
    rl = dble(rlaunch)*h1
    do i = 1,3
       box1(i) = dble(xyzlow(i))*h1
       box2(i) = dble(xyzhih(i))*h1
    end do
    bk = tol_given
    eps = dble(tol)*h1
    bz = zeno_done

    if (zeno_done) then
       capd(1) = dble(cap)*h1
       capd(2) = dble(delta_cap)*h1
       do i = 1,3
          do j = 1,3
             alpha_paper(i,j) = alpha_bongo(i,j) * 4.0 * pi
             delta_paper(i,j) = delta_bongo(i,j) * 4.0 * pi
          end do
       end do

       a11(1) = dble(alpha_paper(1,1)) *h3
       a12(1) = dble(alpha_paper(1,2)) *h3
       a13(1) = dble(alpha_paper(1,3)) *h3
       a21(1) = dble(alpha_paper(2,1)) *h3
       a22(1) = dble(alpha_paper(2,2)) *h3
       a23(1) = dble(alpha_paper(2,3)) *h3
       a31(1) = dble(alpha_paper(3,1)) *h3
       a32(1) = dble(alpha_paper(3,2)) *h3
       a33(1) = dble(alpha_paper(3,3)) *h3
       a11(2) = dble(delta_paper(1,1)) *h3
       a12(2) = dble(delta_paper(1,2)) *h3
       a13(2) = dble(delta_paper(1,3)) *h3
       a21(2) = dble(delta_paper(2,1)) *h3
       a22(2) = dble(delta_paper(2,2)) *h3
       a23(2) = dble(delta_paper(2,3)) *h3
       a31(2) = dble(delta_paper(3,1)) *h3
       a32(2) = dble(delta_paper(3,2)) *h3
       a33(2) = dble(delta_paper(3,3)) *h3
       do i = 1,3
          do j = 1,3
             vae(i,j) = h1*vae(i,j)
             wae(i,j) = h1*wae(i,j)
          end do
       end do
       do i = 1,3
          eigens(i) = eigens(i)*h3
       end do
    end if

    bs = kirk_done
    if (kirk_done) then
       rk(1) = dble(kirk)*h1
       rk(2) = dble(delta_kirk)*h1
       surf(1) = dble(surface)*h2
       surf(2) = dble(delta_surface)*h2
       rg2surfd(1) = dble(rg2surf)*h2
       rg2surfd(2) = dble(delta_rg2surf)*h2
    end if

    bi = rg_done
    if (rg_done) then

       v(1) = dble(volume)*h3
       v(2) = dble(delta_volume)*h3
       rg2intd(1) = dble(rg2int)*h2
       rg2intd(2) = dble(delta_rg2int)*h2
       do i = 1,3
          do i1 = 1,3
             tten(i,i1) = tten(i,i1)*h2
          end do
       end do

       do i = 0,81
          q(i) = q(i)/sngl(h1)
          stfac(i,1) = dble(sq(i))
          stfac(i,2) = dble(dsq(i))
       end do

       !> Find first point at which s(q) drops below 0.5

       l_half_flg = .false.
       do i = 1,81
          if (stfac(i,1) < 0.5d0) then
             ka = i - 1
             kb = i
             l_half_flg = .true.
             exit
          end if
       end do

       if (.not. l_half_flg) then
          write(nzno,"('Unable to bracket L-half')")
          stop 'Unable to bracket L-half'
       end if

       qa(1) = dble(q(ka))
       qa(2) = 0.0d0
       qb(1) = dble(q(kb))
       qb(2) = 0.0d0
       sqa(1) = stfac(ka,1)
       sqa(2) = stfac(ka,2)
       sqb(1) = stfac(kb,1)
       sqb(2) = stfac(kb,2)
       call lighthalf(qa, qb, sqa, sqb, lhalf)

    end if

    bp = span_done
    if (bp) then
       rp(1) = dble(span/2.0)*h1
       rp(2) = dble(delta_span/2.0)*h1
    end if

    bpp = shadow_done
    if (bpp) then
       rap(1) = dble(shadow)*h2
       rap(2) = dble(delta_shadow)*h2
    end if

    !! Record results of all calculations on output file

    call poplar(mz, ms, mi, mp, bl, bk, bs, bi, bt, bm, bw, bc, bbf, &
         &      bz, bp, bpp, rl, eps, capd, a11, a12, a13, a21, a22,&
         &      a23, a31, a32, a33, tten, rk, surf, rg2surfd, v,&
         &      rg2intd, rp, rap, temp, mass, visc, buoy, lunit,&
         &      tunit, munit, vunit, qunit, tae, uae, vae, wae,&
         &      q2pade, eigens, xx, q, stfac, lhalf, box1, box2)

  end subroutine report

  !! ----------------------------------------------------------------

  !> Gets time stamp at end of job and writes it to relevant file.
  !!
  !! @warning makes too many assumptions about string sizes.

  subroutine alldone(endend)

    use zeno_files_data

    !! Declare arguments
    character(len=28) :: endend

    !! Local variables
    character(len=32) :: com
    character(len=25) :: efl

    !> get timestamp at end of job

    com(1:7) = 'date > '
    !!          1234567
    com(8:32) = fefl
    call system(com)
    efl = fefl
    call gettime(efl,endend)

  end subroutine alldone

  !! ----------------------------------------------------------------

  !> Read current time from given filename
  !! 
  !! @param dfl   (in) `string`, input file name
  !! @param start (in) `string`, Unix-style date and trime string.
  !! 
  !! @warning makes assumptions about string lengths
  !!
  !! @warning too closely coupled to the way time and date are processed.

  subroutine gettime(dfl,start)

    use zeno_files_data

    !! Declare arguments
    character(len=25) :: dfl
    character(len=28) :: start    !> Unix-style date and time string.

    open(unit=ndfl, file=dfl, status='old')
    read(ndfl, '(a28)') start
    close(ndfl)

  end subroutine gettime

  !! ----------------------------------------------------------------

  !> Routine assembles one line of output (80 chars).  Where
  !! is the line used?
  !! 
  !! @param id     (in)  ZENO run name, `string`
  !! @param putout (in)  String to insert into `round`
  !! @param round  (out) Output line to print later, `string`
  !! 
  !! @warning uses fixed length strings
  !!
  !! @warning should be renamed to reflect functionality!

  subroutine assemble(id, putout, round)

    !! argument declarations

    character(len=25) :: id             !> ZENO run name?
    character(len=30) :: putout         !> String to insert into `round`?
    character(len=80) :: round          !> Output line to print later?

    !! Local variables
    integer :: i, jix, nper
    character(len=80) :: copy

    do i = 1,80
       round(i:i) = ' '
    end do
    round(1:25) = id
    round(26:26) = '|'
    round(27:56) = putout
    round(57:57) = '|'

    do i = 1,57
       jix = index(round,' ')
       copy(1:jix-1) = round(1:jix-1)
       copy(jix:79) = round(jix+1:80)
       copy(80:80) = ' '
       round = copy
    end do

    nper = index(round,' ') - 1
    do i = nper+1,80
       round(i:i) = round(i-nper:i-nper)
    end do

  end subroutine assemble

  !! ----------------------------------------------------------------

  !> Adds @c + or @c - charactera to `pom` parameter
  !! 
  !! @param kk  (in)  integer vector(3)
  !! @param pom (out) output string filled with '+' and '-'
  !! 
  !! @warning name may change to better reflect its functionality.

  subroutine plus_or_minus(kk, pom)

    !! Declare arguments
    integer, dimension(3) :: kk
    character(len=3)      :: pom

    !! Local variables
    integer :: i

    do i = 1,3
       if (kk(i) == -1) then
          pom(i:i) = '-'
       else if (kk(i) == +1) then
          pom(i:i) = '+'
       else
          stop 'pom pom error'
       end if
    end do

  end subroutine plus_or_minus

  !! ================================================================
  !!
  !! Private routines
  !!
  !! ================================================================

  !> Sends error list to report file (`nzno` I/O unit)
  !! 
  !! @warning name may change to better reflect its functionality.

  subroutine show_errors

    use zeno_warnings, only : ncube, ferr, nell,  rerr
    use zeno_files_data

    if (ncube == 1) then
       write(nzno,"('  Point found inside cube.  Only a cause')")
       write(nzno,"('  for concern if the following number is')")
       write(nzno,"('  not close to 1:', g20.8)") ferr
    end if

    if (nell == 1) then
       write(nzno,"('  Overstretched ellipsoid.  Only a cause')")
       write(nzno,"('  for concern if the following number is')")
       write(nzno,"('  not close to 1:',g20.8)") rerr
    end if

  end subroutine show_errors

  !! ----------------------------------------------------------------

  !> Record results of all calculations on output file.
  !! 
  !! @warning far too many parameters!
  !! @warning name may change to better reflect its functionality.

  subroutine poplar(mz, ms, mi, mp, &
       &            bl, bk, bs, bi, bt, bm, bw, bc, bbf, bz, bp, bpp, &
       &            rl, eps, cap, &
       &            a11, a12, a13, a21, a22, a23, a31, a32, a33, &
       &            tten, rk, surf, rg2surf, v, rg2int, rp, rap, &
       &            temp, mass, eta, buoy, &
       &            lunit, tunit, munit, vunit, qunit, &
       &            tae, uae, vae, wae, q2pade, eigens, &
       &            xx, q, stfac, lhalf, box1, box2)

    !!GCC$ ATTRIBUTES unused :: mz, ms, mi, mp

    use zeno_files_data

    integer :: mz
    integer :: ms
    integer :: mi
    integer :: mp
    logical :: bu
    logical :: bl
    logical :: bk
    logical :: bs
    logical :: bi
    logical :: bt
    logical :: bm
    logical :: bw
    logical :: bc
    logical :: bbf
    logical :: bz
    logical :: bp
    logical :: bpp
    logical :: bv
    real(dp_k) :: rl
    real(dp_k) :: eps
    real(dp_k), dimension(2) :: cap
    real(dp_k), dimension(2) :: a11
    real(dp_k), dimension(2) :: a12
    real(dp_k), dimension(2) :: a13
    real(dp_k), dimension(2) :: a21
    real(dp_k), dimension(2) :: a22
    real(dp_k), dimension(2) :: a23
    real(dp_k), dimension(2) :: a31
    real(dp_k), dimension(2) :: a32
    real(dp_k), dimension(2) :: a33
    real(dp_k), dimension(3,3) :: tten
    real(dp_k), dimension(2) :: rk
    real(dp_k), dimension(2) :: surf
    real(dp_k), dimension(2) :: rg2surf
    real(dp_k), dimension(2) :: v
    real(dp_k), dimension(2) :: rg2int
    real(dp_k), dimension(2) :: rp
    real(dp_k), dimension(2) :: rap
    real(dp_k), dimension(2) :: temp
    real(dp_k), dimension(2) :: mass
    real(dp_k), dimension(2) :: eta
    real(dp_k), dimension(2) :: buoy
    character(len=6) :: lunit
    character(len=6) :: tunit
    character(len=6) :: munit
    character(len=6) :: vunit
    character(len=6) :: qunit
    real(dp_k) :: tae
    real(dp_k), dimension(3) :: uae
    real(dp_k), dimension(3,3) :: vae
    real(dp_k), dimension(3,3) :: wae
    real :: q2pade
    real(dp_k), dimension(3) :: eigens
    real(dp_k), dimension(2) :: xx
    real, dimension(0:81)  :: q
    real(dp_k), dimension(0:81,2) :: stfac
    real(dp_k), dimension(2) :: lhalf
    real(dp_k), dimension(3) :: box1
    real(dp_k), dimension(3) :: box2

    !! Local variables

    character(len=6)  :: aunit, dunit, funit, intunit
    character(len=6)  :: nounit, volunit, svedunit
    character(len=30) :: name
    character(len=2)  :: dig2
    character(len=5), dimension(12) :: nod

    integer :: i, i1, i2

    logical, dimension(12) :: got

    real :: rg2ot
    real(dp_k) :: aa_pol_num, aa_pol_den, aa_pol
    real(dp_k) :: aa_rg2_num, aa_rg2_den, aa_rg2
    real(dp_k) :: ax2, ax1
    real(dp_k) :: exponent
    real(dp_k), dimension(3) :: eigent
    real(dp_k), dimension(2) :: rv, rat
    real(dp_k), dimension(2) :: fcap, fconv, fpercm
    real(dp_k), dimension(2) :: rgsurf, rgint, t2, c0, d
    real(dp_k), dimension(2) :: fric
    real(dp_k), dimension(2) :: q1, pi, six, three
    real(dp_k), dimension(2) :: q2low, q2hih, q2best, sed
    real(dp_k), dimension(2) :: two, five, ol
    real(dp_k), dimension(2) :: rh, t1, trace, four, vh, malf
    real(dp_k), dimension(2) :: sig11, sig22, sig33, signorm
    real(dp_k), dimension(2) :: sigma, etav, ssig, etam, boltz
    real(dp_k), dimension(2) :: pi12
    real(dp_k), dimension(2) :: celsius
    real(dp_k), dimension(2) :: sqsurf      !> Square root of surface are
    real(dp_k), dimension(2) :: cruss, cray !> Russell and Rayleigh aprox. to C
    real(dp_k), dimension(2) :: pre1, pre2
    real(dp_k), dimension(12,2) :: rod

    !! HACK to get around compiler warnings of unused arguments
    integer :: mi_2, mp_2, ms_2, mz_2

    if (.false.) then
       mi_2 = mi
       mp_2 = mp
       ms_2 = ms
       mz_2 = mz
    end if

    !> List of logical variables that control reporting of results
    !!
    !! - `bt` -- temperature was set
    !! - `bm` -- mass was set
    !! - `bv` -- solvent viscosity was set
    !! - `bbf` -- buoyancy factor was set
    !! - `bl` -- launch radius was determined
    !! - `bk` -- skin thickness was determined
    !! - `bz` -- zeno integration succeeded
    !! - `bs` -- surface integration succeeded
    !! - `bi` -- interior integration succeeded
    !! - `bp` -- projection-onto-line integration succeeded
    !! - `bpp` -- projection-onto-plane integration succeeded
    !! - `bu` -- specific, rather than generic, length units are in use

    bu = .not. (lunit == 'L     ')

    call loadconstants(q1, q2low, q2hih, q2pade, q2best, &
         &             pi, two,three, four, five, six, &
         &             boltz, celsius, pre1, pre2, pi12, fpercm)

    call derive(lunit, aunit, volunit, nounit)

    !> must convert temp to kelvin, eta to centipoise
    call parameters(bt, temp, bc, eta, tunit, vunit, celsius)

    open(unit=nznr, file=fznr, status='unknown')

    write(nzno, "('RESULTS OF CALCULATION:')")
    write(nznr, "('RESULTS OF CALCULATION:')")

    if (bl) then
       name = 'launch radius@'
       call writeno(rl, name, lunit)
       name = 'XYZ_low(x)@'
       call writeno(box1(1), name, lunit)
       name = 'XYZ_low(y)@'
       call writeno(box1(2), name, lunit)
       name = 'XYZ_low(z)@'
       call writeno(box1(3), name, lunit)
       name = 'XYZ_hih(x)@'
       call writeno(box2(1), name, lunit)
       name = 'XYZ_hih(y)@'
       call writeno(box2(2), name, lunit)
       name = 'XYZ_hih(z)@'
       call writeno(box2(3), name, lunit)
    end if

    if (bk) then
       name = 'skin thickness@'
       call writeno(eps, name, lunit)
    end if

    if (bt) then
       name = 'temperature@'
       call writeyes(temp, name, tunit)
    end if

    if (bm) then
       name = 'mass@'
       call writeyes(mass, name, munit)
    end if

    bv = .false.
    if (bc) then
       name = 'solvent viscosity (supplied)@'
       call writeyes(eta, name, vunit)
       bv = .true.
    end if

    if ( (.not. bc) .and. bt .and. bw) then
       call compvisc(temp, eta, vunit)
       bv = .true.
       name = 'solvent viscosity (computed)@'
       call writeyes(eta, name, vunit)
    end if

    if (bbf) then
       name = 'buoyancy factor@'
       call writeyes(buoy, name, nounit)
    end if

    if (bz) then
       name = 'capacitance, C@'
       call writeyes(cap, name, lunit)
       if (bu) then
          funit = 'farad '
          if (lunit == 'm     ') then
             fconv(1) = fpercm(1)*1.0d2
             fconv(2) = fpercm(2)*1.0d2
          else if (lunit == 'cm    ') then
             fconv(1) = fpercm(1)
             fconv(2) = fpercm(2)
          else if (lunit == 'nm    ') then
             fconv(1) = fpercm(1)*1.0d-7
             fconv(2) = fpercm(2)*1.0d-7
          else if (lunit == 'A     ') then
             fconv(1) = fpercm(1)*1.0d-8
             fconv(2) = fpercm(2)*1.0d-8
          end if
          call multiply(cap, fconv, fcap)
          name = 'practical capacitance, C@'
          call writeyes(fcap, name, funit)
       end if
       name = 'polarizability 11@'
       call writeyes(a11, name, volunit)
       name = 'polarizability 12@'
       call writeyes(a12, name, volunit)
       name = 'polarizability 13@'
       call writeyes(a13, name, volunit)
       name = 'polarizability 21@'
       call writeyes(a21, name, volunit)
       name = 'polarizability 22@'
       call writeyes(a22, name, volunit)
       name = 'polarizability 23@'
       call writeyes(a23, name, volunit)
       name = 'polarizability 31@'
       call writeyes(a31, name, volunit)
       name = 'polarizability 32@'
       call writeyes(a32, name, volunit)
       name = 'polarizability 33@'
       call writeyes(a33, name, volunit)
       name = 'pol.eigenvalue 1@'
       eigens(1) = eigens(1)*4.0d0*pi(1)
       eigens(2) = eigens(2)*4.0d0*pi(1)
       eigens(3) = eigens(3)*4.0d0*pi(1)
       call writeno(eigens(1), name, volunit)
       name = 'pol.eigenvalue 2@'
       call writeno(eigens(2), name, volunit)
       name = 'pol.eigenvalue 3@'
       call writeno(eigens(3), name, volunit)
       name = 'pol.shape ratio 1@'
       call writeno(xx(1), name, nounit)
       name = 'pol.shape ratio 2@'
       call writeno(xx(2), name, nounit)
       call rudnick_gaspari(eigens, aa_pol_num, aa_pol_den, aa_pol)
       name = 'A(RG).pol.num@'
       call writeno(aa_pol_num, name, aunit)
       name = 'A(RG).pol.den@'
       call writeno(aa_pol_den, name, aunit)
       name = 'A(RG).pol@'
       call writeno(aa_pol, name, nounit)
       call multiply(q1, cap, rh)
       name = 'Rh@'
       call writeyes(rh, name, lunit)
       call summer(a11, a22, t1)
       call summer(t1, a33, trace)
       name = 'Tr(alpha)@'
       call writeyes(trace, name, volunit)
       call divide(trace, three, malf)
       name = '<alpha>@'
       call writeyes(malf, name, volunit)
       name = 'q2(lower)@'
       call writeno(q2low(1), name, nounit)
       name = 'q2(pade)@'
       call writeyes(q2best, name, nounit)
       name = 'q2(pade-no-round)@'
       call writeno(q2best(1), name, nounit)
       name = 'q2(upper)@'
       call writeno(q2hih(1), name, nounit)
       call multiply(q2best, trace, t1)
       call multiply(t1, two, t2)
       call divide(t2, three, t1)
       call divide(t1, five, vh)
       name = 'Vh@'
       call writeyes(vh, name, volunit)
       call multiply(three, vh, t1)
       call divide(t1, four, t2)
       call divide(t2, pi, t1)
       exponent = 1.0d0/3.0d0
       call power(t1, rv, exponent)
       name = 'Rv (viscosity radius)@'
       call writeyes(rv, name, lunit)
       call writehandy(tae, uae, vae, wae, nounit, volunit)
    end if

    if (bp) then
       name = 'LG@'
       call writeyes(rp, name, lunit)
    end if

    if (bpp) then
       name = 'Ap@'
       call writeyes(rap, name, aunit)
       call divide(rap, pi, t1)
       exponent = 1.0d0/2.0d0
       call power(t1, ol, exponent)
       name = 'LO@'
       call writeyes(ol, name, lunit)
    end if

    if (bs) then
       name = 'RK@'
       call writeyes(rk, name, lunit)
       name = 'surface area@'
       call writeyes(surf, name, aunit)
       exponent = 1.0d0/2.0d0
       call power(rg2surf, rgsurf, exponent)
       name = 'Rg (surface)@'
       call writeyes(rgsurf, name, lunit)
       exponent = 1.0d0/2.0d0
       call power(surf, sqsurf, exponent)
       call multiply(sqsurf, pre1, cruss)
       call multiply(sqsurf, pre2, cray)
       name = 'R(Russell)@'
       call writeyes(cruss, name, lunit)
       name = 'R(Rayleigh)@'
       call writeyes(cray, name, lunit)
    end if

    if (bi) then
       name = 'volume@'
       call writeyes(v, name, volunit)
       exponent = 1.0d0/2.0d0
       call power(rg2int, rgint, exponent)
       name = 'Rg (interior)@'
       call writeyes(rgint, name, lunit)
       name = 'L-half@'
       call writeyes(lhalf, name, lunit)
       call ttdiag(tten, eigent)

       do i1 = 1,3
          do i2 = 1,3
             write(dig2,"(i1, i1)") i1, i2

             name = 'Rg2 tensor xx@'
             !!      12345678901234
             name(12:13) = dig2
             call writeno(tten(i1, i2), name, aunit)
          end do
       end do

       name = 'Rg2.eigenvalue 1@'
       call writeno(eigent(1), name, aunit)
       name = 'Rg2.eigenvalue 2@'
       call writeno(eigent(2), name, aunit)
       name = 'Rg2.eigenvalue 3@'
       call writeno(eigent(3), name, aunit)
       ax1 = eigent(2)/eigent(1)
       ax2 = eigent(3)/eigent(2)
       name = 'Rg2.shape ratio 1@'
       call writeno(ax1, name, nounit)
       name = 'Rg2.shape ratio 2@'
       call writeno(ax2, name, nounit)
       call rudnick_gaspari(eigent, aa_rg2_num, aa_rg2_den, aa_rg2)
       name = 'A(RG).rg2.num@'
       call writeno(aa_rg2_num, name, aunit)
       name = 'A(RG).rg2.den@'
       call writeno(aa_rg2_den, name, aunit)
       name = 'A(RG).rg2@'
       call writeno(aa_rg2, name, nounit)
       rg2ot = real(eigent(1) + eigent(2) + eigent(3))
       rg2ot = real(rg2int(1)/rg2ot)
       !! print *,'Rg check should be one:  ', rg2ot

    end if

    if (bi) then
       call multiply(three, v, t1)
       call divide(t1, four, t2)
       call divide(t2, pi, t1)
       exponent = 1.0/3.0
       call power(t1, c0, exponent)
       name = 'C0@'
       call writeyes(c0, name, lunit)
    end if

    if (bz .and. bi) then
       call divide(a11, three, t1)
       call divide(t1, v, sig11)
       call divide(a22, three, t1)
       call divide(t1, v, sig22)
       call divide(a33, three, t1)
       call divide(t1, v, sig33)
       call summer(sig11, sig22, signorm)
       name = '[sigma](xy)@'
       call writeyes(signorm, name, nounit)
       name = '[sigma](z)@'
       call writeyes(sig33, name, nounit)
       call summer(signorm, sig33, sigma)
       name = '[sigma]@'
       call writeyes(sigma, name, nounit)
       call multiply(q2low, sigma, etav)
       name = '[eta](V)(lower)@'
       call writeyes(etav, name, nounit)
       call multiply(q2best, sigma, etav)
       name = '[eta](V)(pade)@'
       call writeyes(etav, name, nounit)
       call multiply(q2hih, sigma, etav)
       name = '[eta](V)(upper)@'
       call writeyes(etav, name, nounit)
       !  recalculate with q2best
       !  so best result gets carried forward
       call multiply(q2best, sigma, etav)
    end if

    if (bi .and. bs) then
       call multiply(six, six, t1)
       call multiply(t1, pi, t2)
       exponent = 1.0d0/3.0d0
       call power(t2, t1, exponent)
       call divide(surf, t1, t2)
       exponent = 2.0d0/3.0d0
       call power(v, t1, exponent)
       call divide(t2, t1, ssig)
       name = 'sphericity@'
       call writeyes(ssig, name, nounit)
    end if

    if (bz .and. bm) then
       call divide(vh, mass, t1)
       call multiply(t1, five, t2)
       call divide(t2, two, t1)
       call conventional(t1, etam, lunit, munit, intunit)
       name = '[eta](M)@'
       call writeyes(etam, name, intunit)
    end if

    if (bz .and. bt .and. bv .and. bu) then
       call multiply(boltz, temp, t1)
       call divide(t1, six, t2)
       call divide(t2, pi, t1)
       call divide(t1, eta, t2)
       call divide(t2, rh, t1)
       call diffuse(t1, d, lunit, dunit)
       name = 'D@'
       call writeyes(d, name, dunit)
    end if

    if (bm .and. bbf .and. bz .and. bu .and. bv) then
       call multiply(mass, buoy, t1)
       call divide(t1, six, t2)
       call divide(t2, pi, t1)
       call divide(t1, eta, t2)
       call divide(t2, rh, t1)
       call svedberg(t1, sed, munit, lunit, svedunit)
       name = 's@'
       call writeyes(sed, name, svedunit)
    end if

    if (bu .and. bv .and. bz) then
       call multiply(six, pi, t1)
       call multiply(t1, eta, t2)
       call multiply(t2, rh, t1)
       call friction(t1, fric, lunit, funit)
       name = 'f@'
       call writeyes(fric, name, funit)
    end if

    if (bi) then
       write(nzno,"(50('='))")
       write(nzno,"('STRUCTURE FACTOR, q, S(q): ')")
       write(nznr,"('STRUCTURE FACTOR, q, S(q): ')")
       do i = 1,81
          call writeme(qunit, nounit, q(i), stfac(i,1), stfac(i,2))
       end do
    end if

    write(nzno,"(50('='))")
    write(nzno,"('DIMENSIONLESS RATIOS:')")
    write(nznr,"('DIMENSIONLESS RATIOS:')")

    do i = 1,2
       rod(1,i)  = cap(i)
       rod(2,i)  = rh(i)
       rod(3,i)  = rv(i)
       rod(4,i)  = rgint(i)
       rod(5,i)  = rgsurf(i)
       rod(6,i)  = lhalf(i)
       rod(7,i)  = c0(i)
       rod(8,i)  = rp(i)
       rod(9,i)  = ol(i)
       rod(10,i) = rk(i)
       rod(11,i) = cruss(i)
       rod(12,i) = cray(i)
    end do

    nod(1) = 'C@'
    nod(2) = 'Rh@'
    nod(3) = 'Rv@'
    nod(4) = 'Rg,i@'
    nod(5) = 'Rg,s@'
    nod(6) = 'Lhlf@'
    nod(7) = 'C0@'
    nod(8) = 'LG@'
    nod(9) = 'LO@'
    nod(10) = 'RK@'
    nod(11) = 'RRus@'
    nod(12) = 'RRay@'

    got(1) = bz
    got(2) = bz
    got(3) = bz
    got(4) = bi
    got(5) = bs
    got(6) = bi
    got(7) = bi
    got(8) = bp
    got(9) = bpp
    got(10) = bs
    got(11) = bs
    got(12) = bs

    do i1 = 1,12
       do i2 = 1,12
          if (i2 /= i1) then
             if (got(i1) .and. got(i2)) then
                call bundle(nod(i1),nod(i2),name)
                t1(1) = rod(i1,1)
                t1(2) = rod(i1,2)
                t2(1) = rod(i2,1)
                t2(2) = rod(i2,2)
                call divide(t1,t2,rat)
                call writeyes(rat,name,nounit)
             end if
          end if
       end do
    end do

    write(nzno,"(50('='))")
    close(nznr)

  end subroutine poplar

  !! ----------------------------------------------------------------

  !> To be filled in.
  !! 
  !! @todo Document parameters
  !! 
  !! @warning name may change to better reflect its functionality.

  subroutine writeno(rl, name, lunit)

    use zeno_files_data

    !! Declare arguments
    real(dp_k) :: rl
    character(len=30) :: name
    character(len=6)  :: lunit

    !! Local variables
    character(len=30) :: line

    call dotout(name, line)

    write(nzno,"((a30), (g15.6), (4x), (a6))") line, rl, lunit
    write(nznr,"((a30), (g15.6), (4x), (a6))") line, rl, lunit

  end subroutine writeno

  !! ----------------------------------------------------------------

  !> Fills the string `line` with `name` followed by repetitions of
  !! the filler @c " ."
  !!
  !! @warning assumes a string of length 30!
  !! @warning uses the special character '@@' as a field marker.

  subroutine dotout(name, line)

    !! Declare arguments
    character(len=30) :: name, line

    !! Local variables
    integer :: i, jj

    jj = index(name,'@')
    line = name
    line(jj:jj) = ' '

    do i = jj+1,30
       if (mod(i,2)==1) line(i:i) = ' '
       if (mod(i,2)==0) line(i:i) = '.'
    end do

  end subroutine dotout

  !! ----------------------------------------------------------------

  !> To be filled in!
  !!
  !! @todo Document parametes
  !! 
  !! @warning name may change to better reflect its functionality.

  subroutine writeyes(tx, name, u)

    use zeno_files_data

    !! Declare arguments
    real(dp_k), dimension(2) :: tx
    character(len=30) :: name
    character(len=6)  :: u

    !! Local variables
    character(len=30) :: line, mane
    character(len=49) :: longname
    integer           :: j, jax

    call longjohn(tx, name, u, mane, line, longname)

    write(nzno,'((a30), (1x), (a49))') line, longname

    jax = index(mane,'@')
    if (jax==0) stop '@ error'
    do j = jax,30
       mane(j:j) = ' '
    end do

    write(nznr,'((a30),(2x),(2g20.8),(2x),(a6))') mane, tx, u

  end subroutine writeyes

  !! ----------------------------------------------------------------

  !> To be filled in
  !! 
  !! @todo Document parameters
  !! 
  !! @warning name may change to better reflect its functionality.

  subroutine writehandy(tae, uae, vae, wae, nounit, vunit)

    use zeno_files_data

    !! Declare arguments
    real(dp_k)                 :: tae
    real(dp_k), dimension(3)   :: uae
    real(dp_k), dimension(3,3) :: vae
    real(dp_k), dimension(3,3) :: wae
    character(len=6)         :: nounit
    character(len=6)         :: vunit

    !! Local variables
    integer :: i, j

    write(nzno,"(('     T!   . . . . '),(g20.10),(1x),(a6))") tae, nounit
    write(nznr,"(('     T!   . . . . '),(g20.10),(1x),(a6))") tae, nounit

    do i = 1,3
       write(nzno,"('     U!',i1,'  . . . . ',g20.10,1x,a6)") i, uae(i), nounit
       write(nznr,"('     U!',i1,'  . . . . ',g20.10,1x,a6)") i, uae(i), nounit
    end do

    do i = 1,3
       do j = 1,3
          write(nzno,"('     V!',2i1,' . . . . ',g20.10,1x,a6)") &
               & i, j, vae(i,j), vunit
          write(nznr,"('     V!',2i1,' . . . . ',g20.10,1x,a6)") &
               & i, j, vae(i,j), vunit
       end do
    end do

    do i = 1,3
       do j = 1,3
          write(nzno,"('     W!',2i1,' . . . . ',g20.10,1x,a6)") &
               & i, j, wae(i,j), vunit
          write(nznr,"('     W!',2i1,' . . . . ',g20.10,1x,a6)") &
               & i, j, wae(i,j), vunit
       end do
    end do

  end subroutine writehandy

  !! ----------------------------------------------------------------

  !> To be filled in
  !! 
  !! @todo Document parameters
  !! 
  !! @warning name may change to better reflect its functionality.

  subroutine writeme(qunit, nounit, q, stfac1, stfac2)

    use zeno_files_data

    !! Declare arguments
    character(len=6) :: qunit, nounit
    real             :: q
    real(dp_k)         :: stfac1, stfac2

    !! Local variables

    real(dp_k), dimension(2) :: t
    character(len=30)      :: name, line, mane
    character(len=49)      :: longname
    integer                :: j, jax

    t(1) = stfac1
    t(2) = stfac2

    write(name,"(g10.5,1x,a6,'@')") q, qunit
    call longjohn(t, name, nounit, mane, line, longname)

    write(nzno,"(a30,1x,a49)") line, longname

    jax = index(mane,'@')
    if (jax==0) stop '@ error'
    do j = jax,30
       mane(j:j) = ' '
    end do

    write(nznr,"(a30,2x,2g20.8,2x,a6)") mane,t,nounit

  end subroutine writeme

  !! ----------------------------------------------------------------

  !> Routine fills a string by first counting the number of decimal
  !! digits in a real number and then, among other things, formatting
  !! the real number into the string
  !!
  !! @warning name may change to better reflect its functionality.

  subroutine longjohn(tx, name, u, mane, line, longname)

    !! Declare arguments
    real(dp_k), dimension(2) :: tx
    character(len=30)      :: name
    character(len=6)       :: u
    character(len=30)      :: mane
    character(len=30)      :: line
    character(len=49)      :: longname

    !! Local variables
    real(dp_k), dimension(2) :: t
    character(len=20)      :: string
    character(len=1)       :: byte
    character(len=6)       :: blanks
    character(len=49)      :: copy
    integer                :: jj, jq, k1, k2, look, ndes, nex

    mane = name                 !  Operate on a copy
    t(1) = tx(1)                !  Operate on a copy
    t(2) = tx(2)

    blanks = '      '

    !> @warning compares two pairs of doubles with ==
    if (t(1) == 0.0d0 .and. t(2) == 0.0d0) then

       k1 = 0
       k2 = 0

       !> @warning compares two doubles with ==
    else if (t(2) == 0.0d0) then

       k1 = nint(t(1))
       k2 = 0
       nex = 0

       do while (k1 > 999999 .or. 99999 >= k1)

          if (k1 > 999999) then
             t(1) = t(1)/10.0d0
             nex = nex + 1
          else if (99999 >= k1) then
             t(1) = t(1)*10.0d0
             nex = nex - 1
          else
             stop 'longjohn 1 error'
          end if

          k1 = nint(t(1))
       end do

    else

       nex = 0

       do while (t(2) < 0.95d0 .or. 9.5 <= t(2))

          if (t(2) < 0.95d0) then
             t(1) = t(1) * 10.0d0
             t(2) = t(2) * 10.0d0
             nex = nex - 1
          else if (9.5d0 <= t(2)) then
             t(1) = t(1) / 10.0d0
             t(2) = t(2) / 10.0d0
             nex = nex + 1
          else
             stop 'longjohn error 3'
          end if

       end do

       k1 = nint(t(1))
       k2 = nint(t(2))

    end if

    write(string,"(i10,'.(',i7,')')") k1,k2

    call pack20(string)
    if (string(1:1) == '-') then
       ndes = 3
    else
       ndes = 2
    end if

    look = index(string,'.')
    if (look == 0) stop 'longjohn error 2'
    do while (look /= ndes)
       if (look > ndes) then
          call shiftleft(string)
          nex = nex + 1
       else
          call shiftright(string)
          nex = nex - 1
       end if
       look = index(string,'.')
    end do

    call dotout(name,line)

    if (nex < 0) then
       byte = '-'
       nex = iabs(nex)
    else
       byte = '+'
    end if

    if (nex == 0) then
       write(longname,"(a20, 22x,'@',a6)") string, blanks
    else if (nex >= 10) then
       write(longname,"(a20, 'E',a1,i20,'@',a6)") string, byte, nex, blanks
    else
       write(longname,"(a20, 'E',a1,i20.2,'@',a6)") string, byte, nex, blanks
    end if

    call pack49(longname)

    jq = index(longname,'.')
    if (longname(jq+1:jq+1) == '(') then
       longname(jq:jq) = ' '
       call pack49(longname)
    end if

    jj = index(longname,'@')
    longname(jj:jj) = ' '
    if (longname(1:1) /= '-') then
       copy(2:49) = longname(1:48)
       copy(1:1) = ' '
       longname = copy
    end if
    longname(19:24) = u

  end subroutine longjohn

  !! ----------------------------------------------------------------

  !> To be filled in!
  !! 
  !! @warning what is special about 20?
  !! @warning name may change to better reflect its functionality.
  
  subroutine pack20(longname)

    !! Declare arugment 
    character(len=20) :: longname

    !! Local variables
    character(len=1) :: b1, b2
    integer          :: nwk

    nwk = 1
    do while (nwk < 20)
       if (nwk == 0) nwk = 1
       b1 = longname(nwk:nwk)
       b2 = longname(nwk+1:nwk+1)
       if (b1 == ' ' .and. b2 /= ' ') then
          longname(nwk:nwk) = b2
          longname(nwk+1:nwk+1) = b1
          nwk = nwk - 1
       else
          nwk = nwk + 1
       end if
    end do

  end subroutine pack20

  !! ----------------------------------------------------------------

  !> To be filled in!
  !! 
  !! @warning what is special about 49?
  !! @warning name may change to better reflect its functionality.

  subroutine pack49(longname)

    !! Declare argument
    character(len=49) :: longname

    !! Local variables
    character(len=1) :: b1,b2
    integer          :: nwk

    nwk = 1
    do while (nwk < 49)
       if (nwk == 0) nwk = 1
       b1 = longname(nwk:nwk)
       b2 = longname(nwk+1:nwk+1)
       if (b1 == ' ' .and. b2 /= ' ') then
          longname(nwk:nwk) = b2
          longname(nwk+1:nwk+1) = b1
          nwk = nwk - 1
       else
          nwk = nwk + 1
       end if
    end do

  end subroutine pack49

  !! ----------------------------------------------------------------

  !> to be filled in!
  !! 
  !! @warning assumes a string of a particular length!
  !! @warning name may change to better reflect its functionality.

  subroutine shiftleft(string)

    !! Declare argument
    character(len=20) :: string

    !! Local variables
    character(len=1) :: b1, b2
    integer          :: idx

    idx = index(string,'.')
    if (idx == 0) stop 'shiftleft error 1'
    if (idx == 1) stop 'shiftleft error 2'

    b1 = string(idx-1:idx-1)
    b2 = string(idx:idx)
    string(idx-1:idx-1) = b2
    string(idx:idx) = b1

  end subroutine shiftleft

  !! ----------------------------------------------------------------

  !> to be filled in!
  !! 
  !! @warning explicit assumption about string length!

  subroutine shiftright(string)

    !! Declare argument
    character(len=20) :: string

    !! Local variables
    character(len=1) :: b1, b2
    integer          :: idx

    idx = index(string,'.')
    if (idx == 0) stop 'shiftright error 1'
    if (idx == 20) stop 'shiftright error 2'

    b1 = string(idx+1:idx+1)
    b2 = string(idx:idx)
    string(idx+1:idx+1) = b2
    string(idx:idx) = b1

  end subroutine shiftright

  !! ----------------------------------------------------------------

  !> Fill `name` with `nod1` and `nod2` along marker.
  !!
  !! @warning uses the special character "@@" as a maker.
  !! @warning assumes particular string lengths!

  subroutine bundle(nod1, nod2, name)

    !! Declare arguments
    character(len=5)  :: nod1, nod2
    character(len=30) :: name

    !! Local variables
    integer :: i, kk, lk

    lk = 1
    kk = 1

    do while (nod1(kk:kk) /= '@')
       name(lk:lk) = nod1(kk:kk)
       kk = kk + 1
       lk = lk + 1
    end do

    name(lk:lk) = '/'

    do i = 1,5
       lk = lk + 1
       name(lk:lk) = nod2(i:i)
    end do

  end subroutine bundle

  !! ================================================================

end module zeno_output

!!! ================================================================

!!! Local Variables:
!!! mode: f90
!!! time-stamp-line-limit: 30
!!! fill-column: 72
!!! End:
