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
!! Main program for ZENO in Fortran 2008!
!! 
!! @author Mark Mansfield (formerly of Stevens Institute) and students
!! @author Walid Keyrouz (walid.keyrouz@nist.gov)
!! @date   Fri Nov  1 13:53:09 2013 EDT
!! 
!! @defgroup main Main Program
!! 
!! Main program for ZENO in Fortran 2008!
!! 
!! Program allocates arrays for the rest of the modules, reads and
!! parses the input file, invokes the relevant computation, and then
!! generates the relevant output report.
!! 
!! @warning Defines too many global variables.
!!
! Time-stamp: <2015-01-05 16:38:26 walid>
! 
! ================================================================

program zeno

  !! ================================================================

  use numeric_kinds

  use zeno_setup, only : setup
  use zeno_parser, only : parse
  use zeno_output, only : begin_output, assemble, alldone, report
  use zeno_problem_geometry, only : do_launch
  use zeno_tensors, only : pade
  use zeno_integrations, only : all_around, wagaroo, blizzard, captain, &
       &                        calipers

  use zeno_debug

  use zeno_warnings, only : ncube, nell
  use zeno_files_data

  use, intrinsic :: iso_fortran_env, only : OUTPUT_UNIT

  implicit none

  integer, parameter :: maxelts = 65000     !> Array capacity
  integer, parameter :: maxwords = 100      !> Maximum number of keywords?

  integer, dimension(maxelts)    :: eltype    !> Elements array
  real, dimension(maxelts,12)    :: bv        !> Elt body values
  real, dimension(maxelts,3,3)   :: rotations !> Elt rotation matrices
  real, dimension(maxelts)       :: saar      !> TBD
  real, dimension(maxelts)       :: vaar      !> TBD
  real(dp_k), dimension(maxelts) :: strikes   !> TBD
  integer                        :: nelts     !> Array size
  character(len=25)              :: id        !> Case ID
  character(len=80)              :: round     !> TBD
  integer, dimension(4)          :: m1        !> TBD

  integer :: mz,mi,ms,mp,m1do                 !> TBD
  real :: tol
  real, dimension(3) :: xyzlow, xyzhih
  real :: cap, delta_cap
  real, dimension(3,3) :: alpha_bongo, delta_bongo
  real :: volume,  delta_volume
  real :: surface, delta_surface
  real :: rg2int,  delta_rg2int
  real(dp_k), dimension(3,3) :: tten
  real, dimension(0:81) :: q, sq, dsq
  real :: rg2surf, delta_rg2surf
  real :: kirk, delta_kirk
  real :: span, delta_span
  logical :: zeno_done, kirk_done, rg_done, span_done
  logical :: launch_done, tol_given, shadow_done
  real(dp_k) :: tae
  real(dp_k), dimension(3) :: uae
  real(dp_k), dimension(3,3) :: vae, wae
  real(dp_k), dimension(3) :: eigens
  real(dp_k), dimension(2) :: xx
  character(len=4) :: actions
  logical :: savehits_z, savehits_i, savehits_s, savehits_p
  character(len=30), dimension(4) :: actout
  character(len=30) :: putout
  character(len=10), dimension(maxwords) :: dictionary
  integer, dimension(maxwords) :: map
  integer :: nwords

  integer :: i, j, nerr
  real    :: delta_shadow, q2pade, rlaunch, rlaunch1, shadow

  !> Code for length units, five allowed values:
  !> 'L ', 'm ', 'cm', 'nm', 'A '
  character(len=2) :: unitcode

  logical :: bt, bm, bw, bc, bbf
  real(dp_k), dimension(2) :: temp, mass, visc !> Input temp, mass, viscosity
  real(dp_k), dimension(2) :: buoy        !> Input buoyancy factor
  character(len=6) :: tunit, munit, vunit !> Units for temp, mass, visc.
  real :: hscale       !> Length scale parameter set by hunits command

  character(len=28) :: start, endend !> Time stamps

  !> The data structures needed for the distance computation are:
  !! - `nneb`, `ninn`
  !! - `nebtab`
  !! - `bubble`
  !! - `bubble_rad`
  !! - `rlist`

  integer :: nneb, ninn
  integer, dimension(maxelts) :: nebtab
  real, dimension(3) :: bubble
  real :: bubble_rad
  real, dimension(maxelts) :: rlist

  !> @verbatim
  !! actions :  Codes for the actions to take:
  !!         :  z = zeno calculation, 
  !!         :  s = surface calculation,
  !!         :  i = interior calculation,
  !!         :  c = C1 interior calculation,
  !!         :  p = project-onto-line calculation
  !! 
  !!         :  Upper-case versions indicate the same thing
  !!         :  but with recording of hits
  !! 
  !!    eltype  :  A code for each element type. 
  !!            :  It maps to one of the nine variables: sphere_code,
  !!            :     triangle_code, disk_code, open_cylinder_code,
  !!            :     solid_cylinder_code, donut_code, ellipsoid_code,
  !!            :     cube_code, pillar_code
  !! 
  !!    nelts   :  The total number of body elements
  !! 
  !!    bv      :  "Body values,"  bv(i,j),j=1,12 are all the
  !!            :  floating point numbers needed to specify the
  !!            :  body
  !! 
  !!    id      :  character string that identifies the body
  !!            :  it is also used to make up the input and output
  !!            :  file names
  !! 
  !!    start   :  time-stamp at beginning of run
  !!    endend  :  time-stamp at end of run
  !! 
  !!    nin     :  file unit for the input/body file
  !! 
  !!    m1(i)   :  Total number of Monte carlo steps to be 
  !!            :  used in the i-th computation.
  !! 
  !!    rotations  :  For disks, cylinders, and tori:
  !!               :  This stores
  !!               :  a rotation matrix required to rotate the
  !!               :  axis of the element into the z-axis.
  !!               :  For ellipsoids, it stores a rotation matrix
  !!               :     required to rotate n1 to x-axis, n2 to 
  !!               :     y-axis, n3 to z-axis
  !! 
  !!    saar(j)   :  The surface area of the j-th element
  !! 
  !!    vaar(j)   :  The volume of the j-th element -- only computed
  !!              :    if subroutine wagaroo is entered
  !! 
  !!    tol   :  Skin thickness
  !! 
  !!    xyzlow, xyzhih    :  These define the coordinates of a box
  !!                      :  that envelops the body.  Used to generate
  !!                      :   interior points.
  !! 
  !!    cap,delta_cap   :  The capacitance and its uncertainty
  !! 
  !!    alpha_bongo,delta_bongo
  !!                     :  The polarizability tensor, and its
  !!                     :  uncertainty
  !! 
  !!    tae,uae,vae,wae are being calculated and reported for
  !!    the cases in which an ensemble average will be taken
  !! 
  !!    eigens          :  The three eigenvalues of the polarizability
  !!    xx              :  The shape-space coordinates of 
  !!                                the polarizability
  !!    q2pade          :  The proportionality between polarizability
  !!                                and intrinsic viscosity, as 
  !!                                determined by Pade approximation
  !!                                from the eigenvalues of the polariz.
  !! 
  !!    volume,delta_volume
  !!                   :  The volume and its uncertainty
  !! 
  !!    surface,delta_surface
  !!                   :  The surface area and its uncertainty
  !! 
  !!    rg2int, delta_rg2int
  !!                   :  The square radius of gyration of the
  !!                   :  interior, and its uncertainty
  !! 
  !!    tten            :  Rg**2 tensor
  !! 
  !!    q               :  A set of momentum-values at which the
  !!                    :   structure factor is calculated
  !! 
  !!    sq, dsq         :  Values of the structure factor and its
  !!                       uncertainty.
  !! 
  !!    rg2surf, delta_rg2surf
  !!                   :  The square radius of gyration of the
  !!                   :  surface, and its uncertainty
  !! 
  !!    kirk, delta_kirk
  !!                   :  The kirkwood estimate of the hydrodynamic
  !!                   :  radius
  !! 
  !!    span, delta_span
  !!                   :  The Giddings length and its uncertainty
  !! 
  !!    zeno_done       :  The zeno calculation terminated successfully
  !! 
  !!    kirk_done       :  The surface calculation terminated successfully
  !! 
  !!    rg_done         :  The interior calculation terminated successfully
  !! 
  !!    span_done       :  The project-onto-line 
  !!                            calculation terminated successfully
  !! 
  !!    shadow_done     :  The projection-onto-plane calculation was done
  !! 
  !!    launch_done     :  The launch sphere has been generated
  !! 
  !!    tol_given       :  The skin thickness was supplied in the bod file
  !! 
  !!    ncube,ferr      :  Sometimes, when using the cube body element,
  !!                    :  some points can be found "slightly" inside the
  !!                    :  cube.  If these values get set -- the user
  !!                    :  will be warned by an error that this has happened,
  !!                    :  and the user will be told how far inside the
  !!                    :  cube the offending point(s) were found.
  !! 
  !!    nell,rerr       :  I have put in a trap to check for overstretching
  !!                   :  of ellipsoids, because I was not completely
  !!                   :  confident of the stretching equations.  If these
  !!                   :  values get set it is an indication that over-
  !!                   :  stretching has happened.  The user will be
  !!                   :  warned about this with an error statement.
  !! 
  !!    mz, mi, ms, mp  :  Actual computation length of each integration
  !! 
  !!    bt      : temperature was set in body file
  !!    bm      : mass was set in the body file
  !!    bw      : solvent=water command was found in the body file
  !!    bc      : viscosity supplied in body file
  !!    bbf     : buoyancy factor was set in the body file
  !! 
  !!    dictionary:     :  List of all the recognized words in the body file
  !!    map:            :  Integer code corresponding to each entry in
  !!                            dictionary
  !!    nwords:         :  Length of dictionary = total number of words
  !!                            that are valid entries in the body file
  !! @endverbatim

  !> Execution starts here

  ncube = 0
  nell = 0
  nerr = 0

  launch_done = .false.
  zeno_done = .false.
  kirk_done = .false.
  rg_done = .false.
  span_done = .false.
  shadow_done = .false.

  !> Parse the invocation string:
  call setup(id, m1, actions, start, actout, &
       &     savehits_z, savehits_i, savehits_s, savehits_p)

  !> Parse the body file:
  call parse(maxelts, nelts, eltype, bv, tol, rotations, tol_given, &
       &     unitcode, bt, bm, bw, bc, bbf, temp, tunit, mass, munit, &
       &     visc, vunit, buoy, hscale, rlaunch, launch_done, &
       &     dictionary, map, maxwords, nwords)

  !> Begin output file:
  call begin_output(id, nelts, start)

  !> Work out the radius of the launch sphere and the enveloping box
  call do_launch(maxelts, nelts, eltype, bv, rlaunch1, &
       &         rotations, xyzlow, xyzhih)

  !> 
  !! `rlaunch` = user-supplied launch radius
  !! 
  !! `rlaunch1` = launch radius computed from the structure
  !! 
  !! At this point `launch_done = .true.` if the user supplied the
  !! launch radius in the body file.
  !! 
  !! The user-supplied launch radius is used if it is equal to or
  !! greater than the computed radius.  An error exit occurs if the
  !! user-supplied radius is less than the computed radius.
  !! 
  !! If no launch radius was supplied, we use the computed radius.

  if (launch_done) then
     if (rlaunch1 > rlaunch) then
        write(nzno,"('  Invalid launch radius specified.')")
        stop 'Invalid launch radius'
     end if
  else
     rlaunch = rlaunch1
     launch_done = .true.
  end if

  !> At this point, rlaunch contains the launch radius which will be
  !! used, and launch_done = .true. always.

  if (launch_done) then
     if (.not. silent) then
        write(*,"('Launch radius = ',(g15.7))") rlaunch
        write(*,"('XYZ(low)  = ',3(g15.7))") (xyzlow(j),j=1,3)
        write(*,"('XYZ(high) = ',3(g15.7))") (xyzhih(j),j=1,3)
     end if
  end if

  do i = 1,4
     m1do = m1(i)
     putout = actout(i)
     call assemble(id,putout,round)

     if (actions(i:i) == 'i' .or. actions(i:i) == 'I') then

        !> Do the interior integrations:

        call all_around(maxelts, eltype, bv, nelts, &
             &          xyzlow, xyzhih, rlaunch, rotations, m1do, &
             &          rg2int, delta_rg2int, volume, delta_volume, &
             &          rg_done, mi, id, q, sq, dsq, round, savehits_i, tten)

     else if (actions(i:i) == 'c' .or. actions(i:i) == 'C') then

        !> Do the covered-interior integrations:

        call wagaroo(maxelts, eltype, bv, nelts, rlaunch, rotations, m1do, &
             &       rg2int, delta_rg2int, volume, delta_volume, rg_done, &
             &       mi, id, q, sq, dsq, round, savehits_i, vaar, tten)

     else if (actions(i:i) == 'z' .or. actions(i:i) == 'Z') then

        !> Do the zeno integrations:

        !> 
        !! This value of nneb optimizes the run time for about 10
        !! different random coil models

        nneb = 3

        call blizzard(maxelts, eltype, bv, nelts, m1do, tol, rlaunch, &
             &        rotations, cap, alpha_bongo, tol_given, zeno_done, &
             &        delta_cap, delta_bongo, mz, id, tae, uae, vae, wae, &
             &        round, bubble, bubble_rad, nebtab, nneb, ninn, rlist, &
             &        strikes, savehits_z)

        call pade(alpha_bongo, q2pade, eigens, xx)

     else if (actions(i:i) == 's' .or. actions(i:i) == 'S') then

        !> Do the surface integrations:

        call captain(maxelts, eltype, bv, nelts, m1do, rotations, &
             &       kirk_done, saar, kirk, delta_kirk, &
             &       surface, delta_surface, rg2surf, delta_rg2surf, &
             &       ms, id, round, savehits_s)
    
     else if (actions(i:i) == 'p' .or. actions(i:i) == 'P') then
            
        !> Do the project-onto-line integration:

        call calipers(maxelts, eltype, bv, nelts, m1do, rotations, &
             &        span_done, span, delta_span, shadow_done, shadow, &
             &        delta_shadow, rlaunch, mp, id, round, savehits_p)

     end if

  end do

  call alldone(endend)

  call report(id, actions, m1, nelts, tol, rlaunch, &
       &      cap, delta_cap, alpha_bongo, delta_bongo, tten, &
       &      volume, delta_volume, surface, delta_surface, &
       &      rg2int, delta_rg2int, rg2surf, delta_rg2surf, &
       &      kirk, delta_kirk, span, delta_span, &
       &      shadow, delta_shadow, shadow_done, &
       &      zeno_done, kirk_done, rg_done, span_done, &
       &      launch_done, tol_given, unitcode, mz, mi, ms, mp, &
       &      bt, bm, bw, bc, bbf, temp, tunit, mass, &
       &      munit, visc, vunit, buoy, hscale, &
       &      tae, uae, vae, wae, q2pade, eigens, xx, &
       &      q, sq, dsq, xyzlow, xyzhih)

  if (.not. silent) then
     write(nzno,"('END:   ',(a28))") endend
     write(*,   "('END:   ',(a28))") endend
  end if

  close(nzno)

  if (.not. silent) then
     call write_RNG_log(OUTPUT_UNIT)
     call write_LS_log(OUTPUT_UNIT)
     call write_DIST_log(OUTPUT_UNIT)

     call write_RNG_log(nlog)
     call write_LS_log(nlog)
     call write_DIST_log(nlog)
  end if

  close(nlog)

  stop

end program zeno

!!! ================================================================

!!! Local Variables:
!!! mode: f90
!!! time-stamp-line-limit: 30
!!! fill-column: 72
!!! End:
