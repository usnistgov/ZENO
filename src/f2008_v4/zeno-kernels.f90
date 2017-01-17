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
!! This file provides the `zeno_kernels` module which defines routines
!! that perform kernel computations.
!! 
!! @author Mark Mansfield (formerly of Stevens Institute) and students
!! @author Walid Keyrouz (walid.keyrouz@nist.gov)
!! @date   Tue Sep 24 10:18:07 2013 EDT
!!
!! @defgroup zeno_kernels Kernel Computations
!! 
!! Groups computational routines that seem to be ZENO specific.  It
!! defines the following public entry points: `lighthalf`, `compvisc`
!! 
! Time-stamp: <2015-02-23 10:53:54 wtk>
!
! ================================================================

module zeno_kernels

  !! ================================================================

  use numeric_kinds
  use numeric_constants

  use zeno_math_utils

  implicit none
  private
  public lighthalf, compvisc, loadconstants, rudnick_gaspari

contains

  !! ================================================================
  !!
  !! Public routines
  !!
  !! ================================================================

  !> @brief Linear interpolation with error propagation to estimate
  !! L-half

  subroutine lighthalf(qa, qb, sqa, sqb, lhalf)

    !! Declare arguments

    real(dp_k), dimension(2) :: qa        !> TBD
    real(dp_k), dimension(2) :: qb        !> TBD
    real(dp_k), dimension(2) :: sqa       !> TBD
    real(dp_k), dimension(2) :: sqb       !> TBD
    real(dp_k), dimension(2) :: lhalf     !> TBD

    !! Local variables
    real(dp_k), dimension(2) :: mhalf, mqa, msqb
    real(dp_k), dimension(2) :: t1, t2, t3, t4, t5, t6

    mhalf(1) = -0.5d0
    mhalf(2) = 0.0d0
    mqa(1) = -qa(1)
    mqa(2) = qa(2)
    msqb(1) = -sqb(1)
    msqb(2) = sqb(2)

    call summer(qb, mqa, t1)
    call summer(sqa, mhalf, t2)
    call summer(sqa, msqb, t3)
    call multiply(t1, t2, t4)
    call divide(t4, t3, t5)
    call summer(qa, t5, t6)
    t1(1) = 1.0d0
    t1(2) = 0.0d0
    call divide(t1, t6, lhalf)

  end subroutine lighthalf

  !! ----------------------------------------------------------------

  !> @brief Compute viscosity of water as a function of temperature.
  !! 
  !! Compute viscosity of water as a function of temperature.  Return
  !! water viscosity in `eta` and unit used in `vunit`
  !! 
  !! @param temp  (in)  temperature, double precision vector (2)
  !! @param eta   (out) water viscosity, double precision vector (2)
  !! @param vunit (out) unit of viscosity, string(length <= 6)
  !! 
  !! @see CRC handbook of chemistry and physics, 55th edition p. F49
    
  subroutine compvisc(temp, eta, vunit)

    !! Declare arguments
    real(dp_k), dimension(2) :: temp     !> temperature
    real(dp_k), dimension(2) :: eta      !> water viscosity
    character(len=6)         :: vunit    !> unit of viscosity

    !! Local variables
    real(dp_k) :: tc, tm, den, z, num, tp

    !> Compute viscosity of water as a function of temperature
    !! 
    !! Ref. CRC handbook of chemistry and physics, 55th edition p. F49

    tc = temp(1) - 273.15d0

    if (tc <= 20.0d0) then
       tm = tc - 20.0d0
       den = 998.333d0 + 8.1855d0*tm + 0.00585d0*(tm**2)
       z = 1301.0d0/den - 3.30233d0
       eta(1) = 10.0d0**z
       eta(1) = eta(1)*100.0d0
       eta(2) = 0.001d0
    else
       tm = tc - 20.0d0
       tp = tc + 105.0d0
       num = -1.3272d0*tm - 0.001053*tm*tm
       z = num/tp
       eta(1) = 1.002d0 * 10.0d0**z
       eta(2) = 0.0001d0
    end if

    vunit = 'cp    '

  end subroutine compvisc

  !! ----------------------------------------------------------------

  !> Straight-forward code for setting values of parameters with numeric
  !! constants.
  !! 
  !! @warning Far too many parameters!
  !! 
  !! @warning Consider using compile-time constant definitions.

  subroutine loadconstants(q1, q2low, q2hih, q2pade, q2best, &
       &                   pi, two, three, four, five, six, &
       &                   boltz, celsius, pre1, pre2, pi12, fpercm)

    !! Declare arguments
    real(dp_k), dimension(2) :: q1
    real(dp_k), dimension(2) :: q2low
    real(dp_k), dimension(2) :: q2hih
    real                     :: q2pade
    real(dp_k), dimension(2) :: q2best
    real(dp_k), dimension(2) :: pi
    real(dp_k), dimension(2) :: two, three, four, five, six
    real(dp_k), dimension(2) :: boltz
    real(dp_k), dimension(2) :: celsius
    real(dp_k), dimension(2) :: pre1, pre2, pi12
    real(dp_k), dimension(2) :: fpercm

    q1(1) = 1.0d0
    q1(2) = 0.01d0

    q2low(1) = 3.0d0/4.0d0
    q2low(2) = 0.0d0

    q2hih(1) = 5.0d0/6.0d0
    q2hih(2) = 0.0d0

    q2best(1) = q2pade
    q2best(2) = q2best(1)*(0.015d0) 
        
    pi(1) = M_PI
    pi(2) = 0.0d0

    six(1) = 6.0d0
    six(2) = 0.0d0

    three(1) = 3.0d0
    three(2) = 0.0d0
        
    four(1) = 4.0d0
    four(2) = 0.0d0

    two(1) = 2.0d0
    two(2) = 0.0d0

    five(1) = 5.0d0
    five(2) = 0.0d0

    boltz(1) = 1.38065d-23
    boltz(2) = 1.0d-23/1.0d5

    celsius(1) = 273.15d0
    celsius(2) = 0.0d0

    pre1(1) = dsqrt(1.0d0/(4.0d0*pi(1)))
    pre2(1) = dsqrt(1.0d0/pi(1))*(2.0d0/pi(1))
    pre1(2) = 0.0d0
    pre2(2) = 0.0d0

    pi12(1) = (12.0d0*pi(1))**(1.0d0/3.0d0)
    pi12(2) = 0.0d0

    fpercm(1) = 8.987d0
    fpercm(2) = 0.001d0
    fpercm(1) = fpercm(1)*1.0d11
    fpercm(2) = fpercm(2)*1.0d11

  end subroutine loadconstants

  !! ----------------------------------------------------------------

  !> @brief Subroutine `rudnick_gaspari` doc to be filled in.

  subroutine rudnick_gaspari(ev, aanum, aaden, aa)

    !! Declare arguments
    real(dp_k), dimension(3) :: ev
    real(dp_k) :: aanum, aaden, aa

    !! Local variables
    real(dp_k) :: ss, mm

    ss = ev(1) + ev(2) + ev(3)
    mm = ev(1)*ev(2) + ev(1)*ev(3) + ev(2)*ev(3)
    aanum = ss*ss - 3.0d0*mm
    aaden = ss*ss
    aa = aanum/aaden
 
  end subroutine rudnick_gaspari

  !! ================================================================

end module zeno_kernels

!!! ================================================================

!!! Local Variables:
!!! mode: f90
!!! time-stamp-line-limit: 30
!!! fill-column: 72
!!! End:
