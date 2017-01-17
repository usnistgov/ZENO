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
!! This file provides the `zeno_units` module which defines routines that
!! manipulate units.
!! 
!! @author Mark Mansfield (formerly of Stevens Institute) and students
!! @author Walid Keyrouz (walid.keyrouz@nist.gov)
!! @date   Fri Sep 27 15:36:10 2013 EDT
!!
!! @defgroup zeno_units Units
!! 
!! Groups routines that manipulate *units*.  It defines the following
!! entry points: `parameters`, `derive`.
!! 
!! @todo
!! Consider having an explicit representation of `units`.
!! 
! Time-stamp: <2015-01-05 16:45:29 walid>
!
! ================================================================

module zeno_units

  !! ================================================================

  use numeric_kinds
  use zeno_math_utils

  implicit none
  private
  public conventional, derive, diffuse, friction, parameters, svedberg

contains

  !! ================================================================
  !!
  !! Public routines
  !!
  !! ================================================================

  !> Convert units in parameters: temp to kelvin, eta to centipoise
  !! 
  !! @warning needs to be renamed to reflect functionality!

  subroutine parameters(bt, temp, bc, eta, tunit, vunit, celsius)

    !! Declare arguments
    logical                  :: bt
    real(dp_k), dimension(2) :: temp
    logical                  :: bc
    real(dp_k), dimension(2) :: eta
    character(len=6)         :: tunit
    character(len=6)         :: vunit
    real(dp_k), dimension(2) :: celsius

    !! Local variables
    real(dp_k), dimension(2) :: t1, t2

    !  must convert temp to kelvin, eta to centipoise

    if (bt) then
       t1(1) = temp(1)
       t1(2) = temp(2)
       if (tunit == 'C     ') then
          t2(1) = celsius(1)
          t2(2) = celsius(2)
          tunit = 'K     '
       else if (tunit == 'K     ') then
          t2(1) = 0.0d0
          t2(2) = 0.0d0
       end if
       call summer(t1,t2,temp)
    end if

    if (bc) then
       t1(1) = eta(1)
       t1(2) = eta(2)
       if (vunit == 'cp    ') then
          eta(1) = t1(1)
          eta(2) = t1(2)
       else if (vunit == 'p     ') then
          t2(1) = 100.0d0
          t2(2) = 0.0d0
          call multiply(t1,t2,eta)
          vunit = 'cp    '
       end if
    end if

  end subroutine parameters

  !! ----------------------------------------------------------------

  !> Formats "exponent" for length unit?
  !! 
  !! @warning should be renamed to reflect functionality!

  subroutine derive(lunit, aunit, volunit, nounit)

    !! Declare arugments
    character(len=6) :: lunit
    character(len=6) :: aunit, volunit, nounit

    nounit = '      '
    aunit = lunit
    volunit = lunit

    if (aunit(2:2) == ' ') then
       aunit(2:3) = '^2'
       volunit(2:3) = '^3'
    else
       aunit(3:4) = '^2'
       volunit(3:4) = '^3'
    end if

  end subroutine derive

  !! ----------------------------------------------------------------

  !> @brief Converts to cm^3/gm and returns result in `etam`
  !! 
  !! This subroutine does the conversion to cm**3/gm and returns the
  !! result in `etam`.
  !! 
  !! - Possibilities for lunit are  m,cm,nm,A,L
  !! - Possibilities for munit are  kg,g,Da,kDa
  !! 
  !! @pre On entry, `t1` contains the [eta]M value, but in units
  !! lunit**3/munit.

  subroutine conventional(t1, etam, lunit, munit, intunit)

    !! Declare arguments

    real(dp_k), dimension(2) :: t1, etam
    character(len=6)         :: lunit, munit, intunit

    !! Local variables

    real(dp_k), dimension(2) :: t2, conv

    if (lunit == 'L     ') then
       intunit(1:4) = 'L^3/'
       intunit(5:6) = munit(1:2)
       etam(1) = t1(1)
       etam(2) = t1(2)
       return
    end if

    if (lunit == 'm     ') then
       conv(1) = 1.0d6
    else if (lunit == 'cm    ') then
       conv(1) = 1.0d0
    else if (lunit == 'nm    ') then
       conv(1) = 1.0d-21
    else if (lunit == 'A     ') then
       conv(1) = 1.0d-24
    else
       stop 'bad lunit in conventional'
    end if
    conv(2) = 0.0d0

    call multiply(t1, conv, t2)

    if (munit == 'kg    ') then
       conv(1) = 1.0d-3
       conv(2) = 0.0d0
    else if (munit == 'g     ') then
       conv(1) = 1.0d0
       conv(2) = 0.0d0
    else if (munit == 'Da    ') then
       conv(1) = 6.02214d23
       conv(2) = 1.0d18
    else 
       conv(1) = 6.02214d20
       conv(2) = 1.0d15
    end if

    call multiply(t2,conv,etam)

    intunit = 'cm^3/g'

  end subroutine conventional

  !! ----------------------------------------------------------------

  !> Converts the friction constant to (d.s/cm) = dyne.sec/cm and
  !! returns the result in fric
  !! 
  !! @pre
  !! on entry, `t1` contains the friction coefficient in units
  !! (cp lunit)

  subroutine friction(t1, fric, lunit, funit)

    !! Declare arguments
    real(dp_k), dimension(2) :: t1
    real(dp_k), dimension(2) :: fric
    character(len=6)         :: lunit
    character(len=6)         :: funit

    !! Local variable
    real(dp_k), dimension(2) :: conv

    if (lunit == 'm     ') then
       conv(1) = 1.0d0
    else if (lunit == 'cm    ') then
       conv(1) = 1.0d-2
    else if (lunit == 'nm    ') then
       conv(1) = 1.0d-9
    else if (lunit == 'A     ') then
       conv(1) = 1.0d-10
    else
       stop 'bad lunit in friction'
    end if
    conv(2) = 0.0d0

    call multiply(t1, conv, fric)

    funit = 'd.s/cm'

  end subroutine friction

  !! ----------------------------------------------------------------

  !> Converts the sediment coefficient to `svedbergs`.  Upon entry, @c
  !! t1 contains the sedimentation coefficient in units of
  !! (munit)*(centipoise^-1)*(lunit^-1).
  !! 
  !! - Possibilities for `lunit` are m, cm, nm, A
  !! - Possibilities for `munit` are kg, g, Da, kDa
  !! 
  !! @pre
  !! Upon entry, `t1` contains the sedimentation coefficient in units
  !! of (munit)*(centipoise^-1)*(lunit^-1).

  subroutine svedberg(t1, sed, munit, lunit, svedunit)

    !! Declare arguments
    real(dp_k), dimension(2) :: t1        !> value in non-converted units
    real(dp_k), dimension(2) :: sed       !> value in `svedberg` units
    character(len=6)         :: munit
    character(len=6)         :: lunit
    character(len=6)         :: svedunit

    !! Local variables
    real(dp_k), dimension(2) :: conv, t2

    conv(1) = 100.0d0
    conv(2) = 0.0d0
    call multiply(t1,conv,t2)

    !>  Units are now (munit)*(poise^-1)*(lunit^-1) =
    !>  (munit)*(lunit^-1)*(sec)*(cm)*(gm^-1)

    if (munit == 'kg    ') then !  Divide by kg/g
       conv(1) = 1.0d-3
       conv(2) = 0.0d0
    else if (munit == 'g     ') then    !  Divide by g/g
       conv(1) = 1.0d0
       conv(2) = 0.0d0
    else if (munit == 'Da    ') then    !  Divide by Da/g
       conv(1) = 6.02214d23
       conv(2) = 1.0d18
    else if (munit == 'kDa   ') then    !  Divide by kDa/g
       conv(1) = 6.02214d20
       conv(2) = 1.0d15
    else
       print *,'munit = ',munit
       stop 'unrecognized munit in svedberg'
    end if

    call divide(t2,conv,t1)

    !> Units are now (lunit^-1)(sec)(cm)

    if (lunit == 'm     ') then         !  multiply by 1D-2 m/cm
       conv(1) = 1.0d-2
       conv(2) = 0.0d0
    else if (lunit == 'cm    ') then    !  multiply by 1 cm/cm
       conv(1) = 1.0d0
       conv(2) = 0.0d0
    else if (lunit == 'nm    ') then    ! multiply be 1D7 nm/cm
       conv(1) = 1.0d7
       conv(2) = 0.0d0
    else if (lunit == 'A     ') then    ! multiply by 1D8 A/cm
       conv(1) = 1.0d8
       conv(2) = 0.0d0
    else 
       print *,'lunit = ',lunit
       stop 'unrecognized lunit in Svedverg'
    end if

    call multiply(t1,conv,t2)

    !>  Units are now sec, multiply by 1D13 Svedbergs/sec

    conv(1) = 1.0d13
    conv(2) = 0.0d0
    call multiply(t2,conv,sed)
    svedunit = 'Sved  '

  end subroutine svedberg

  !! ----------------------------------------------------------------

  !> @brief converts diffusion constant to (cm^2/s) units.
  !! 
  !! This converts the diffusion constant to (cm^2/s) and returns the
  !! result in `d`.  On entry, t1 contains the diffusion constant in
  !! units J/(cp lunit).
  !! 
  !! @pre On entry, `t1` contains diffusion constant in
  !! "J/(cp lunit) units

  subroutine diffuse(t1, d, lunit, dunit)

    !! Declare arguments
    real(dp_k), dimension(2) :: t1
    real(dp_k), dimension(2) :: d
    character(len=6)         :: lunit
    character(len=6)         :: dunit

    !! Local variable
    real(dp_k), dimension(2) :: conv

    if (lunit == 'm     ') then
       conv(1) = 1.0d7
    else if (lunit == 'cm    ') then
       conv(1) = 1.0d9
    else if (lunit == 'nm    ') then
       conv(1) = 1.0d16
    else if (lunit == 'A     ') then
       conv(1) = 1.0d17
    else
       stop 'bad lunit in diffuse'
    end if
    conv(2) = 0.0d0

    call multiply(t1,conv,d)

    dunit = 'cm^2/s'

  end subroutine diffuse

  !! ================================================================

end module zeno_units

!!! ================================================================

!!! Local Variables:
!!! mode: f90
!!! time-stamp-line-limit: 30
!!! fill-column: 72
!!! End:
