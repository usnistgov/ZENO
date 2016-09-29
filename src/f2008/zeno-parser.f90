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
!! This file defines the `zeno_parser` module which groups functions and
!! routines for parsing ZENO input files.
!!
!! @author Mark Mansfield (formerly of Stevens Institute) and students
!! @author Walid Keyrouz (walid.keyrouz@nist.gov)
!! @date   Wed Sep 18 14:38:05 2013 EDT
!! 
!! @defgroup zeno_parser Parsing for ZENO
!! 
!! Module groups together the subroutines that parse input files.
!! This module provides a single entry point, subroutine `parse`.
!! The other subroutines are all internal routines.  They are:
!! `webster`, `commaparse`, `modifiers`, `nextstring`, 
!! `packspaces`, `cleartabs`, `floatstring`, `pack80`, and 
!! `testword`.
!! 
!! @todo Consider replacing module by C/C++ code that is
!! automatically generated from a grammar definition of the input
!! files.
!! 
! Time-stamp: <2015-01-05 16:44:24 walid>
! 
! ================================================================

module zeno_parser

  !! ================================================================

  use numeric_kinds
  use zeno_vectors
  use zeno_sphere
  use zeno_zeerot

  implicit none
  private

  public parse

contains

  !! ================================================================
  !! 
  !! Public routine
  !! 
  !! ================================================================

  !> @brief Parse the body file, one line at a time.
  !!
  !! @warning too many parameters!  To be replaced by using a parser
  !! generator with a properly defined grammar.
  !!
  !! The list of recognized commands is encoded in an array of
  !! strings, `dictionary`, and a corresponding array of integer
  !! codes, `map`.
  !! 
  !! @verbatim
  !! - SPHERE, sphere, S, s    (sphere_code)
  !!   cx, cy, cz, r
  !!
  !! - TRIANGLE, triangle, T, t  (triangle_code)
  !!   v1x,v1y,v1z,v2x,v2y,v2z,v3x,v3y,v3z  
  !!
  !! - DISK, disk, D, d   (disk_code)
  !!   cx,cy,cz,nx,ny,nz,r 
  !!
  !! - CYLINDER, cylinder    (open_cylinder_code)
  !!   cx,cy,cz,nx,ny,nz,r,l 
  !!
  !! - SOL_CYL, sol_cyl, SC, sc    (solid_cylinder_code)
  !!   cx,cy,cz,nx,ny,nz,r,l 
  !!
  !! - TORUS, torus, TO, to     (donut_code)
  !!   cx,cy,cz,nx,ny,nz,r1,r2 
  !!
  !! - ELLIPSOID, ellipsoid, E, e  (ellipsoid_code)
  !!   cx,cy,cz,n1x,n1y,n1z,n2x,n2y,n2z,aa,bb,cc 
  !!
  !! - CUBE, cube                 (cube_code)
  !!   cx,cy,cz,s 
  !!
  !! - PILLAR, pillar, PI, pi    (pillar_code)
  !!   x1,y1,z1l,z1h,x2,y2,z2l,z2h,x3,y3,z3l,z3h
  !!
  !! - ST, st         (skin_code)
  !!   tol
  !!
  !! - RLAUNCH, rlaunch  (User-defined launch radius, to use
  !!                    in place of radius determined
  !!
  !!  (rlaunch_code)
  !!
  !! - UNITS, units     (units_code)
  !! 
  !!   Modifier      Internal    
  !!   string         code            Meaning
  !!   -------       --------          -------
  !!     m             meter_code     meters
  !!     cm            cm_code        centimeters
  !!     nm            nm_code        nanometers
  !!     A             angstrom_code  Angstroms
  !!     L             length_code    generic or unspecified length units
  !!
  !! - TEMP, temp    (temp_code)
  !!
  !!   First modifier:   number (value of temperature)
  !!
  !!   Second 
  !!   modifier      Internal
  !!   string          code            Meaning
  !!   -------       --------          -------
  !!      C            celcius_code    Celcius
  !!      K            kelvin_code     Kelvin
  !!
  !! - MASS, mass    (mass_code)
  !!
  !!   First modifier:  number (value of mass)
  !!
  !!   Second 
  !!   modifier      Internal
  !!   string          code         Meaning
  !!   -------       --------       -------
  !!      Da           da_code      daltons
  !!      kDa          kda_code     kilodaltons
  !!      g            gram_code    grams
  !!      kg           kg_code      kilograms
  !!
  !! - VISCOSITY, viscosity (visc_code)
  !! 
  !!   Fist modifier:  number (value of solvent viscosity)
  !! 
  !!   Second 
  !!   modifier      Internal
  !!   string          code         Meaning
  !!   -------       --------       -------
  !!       p          poise_code      poise
  !!      cp          cp_code         centipoise
  !!
  !! - BF, bf (bf_code)
  !! 
  !!   First modifier:  number (value of buoyancy factor)
  !! @endverbatim

  subroutine parse(maxelts, nelts, eltype, bv, tol, &
       &           rotations, tolset, unitcode, bt, bm, bw, bc, bbf, &
       &           temp, tunit, mass, munit, visc, vunit, buoy, &
       &           hscale, rlaunch, launch_done, &
       &           dictionary, map, maxwords, nwords)

    use zeno_codes_data
    use zeno_files_data

    !! argument declarations

    integer, intent(in)  :: maxelts           !> capacity of array
    integer, intent(out) :: nelts             !> number of elements processed
    integer, dimension(maxelts)     :: eltype !> type array of all elements
    !> array w/1 row per elt, set from `modifiers` via local array
    !> `download`
    real,    dimension(maxelts, 12) :: bv
    real,    intent(out)   :: tol             !> read, part of ST command?
    real,    dimension(maxelts,3,3) :: rotations
    logical, intent(out)   :: tolset          !> `true` iff `tol` read
    character(len=2)       :: unitcode        !> two char length unit symbol
    logical, intent(out)   :: bt              !> modifier flag for Temp
    logical, intent(out)   :: bm              !> modifier flag for Mass
    logical, intent(out)   :: bw              !> modifier flag for Solvent
    logical, intent(out)   :: bc              !> modifier flag for Viscosity
    logical, intent(out)   :: bbf             !> modifier flag for `bf_code`
    real(dp_k), dimension(2) :: temp          !> temperatures via modifiers
    character(len=6)       :: tunit           !> K or C padded, via modifiers
    real(dp_k), dimension(2) :: mass          !> mass via modifiers
    character(len=6)       :: munit           !> Da, kDa, g, or kg padded,
                                              !> via modifiers
    real(dp_k), dimension(2) :: visc          !> viscosity via modifiers
    character(len=6)       :: vunit           !> p or cp padded, via modifiers
    real(dp_k), dimension(2) :: buoy          !> buoancy (bf) via modifiers
    real, intent(out)      :: hscale          !> 1.0 default, changed by HUNITS
    real, intent(out)      :: rlaunch         !> launch radius, via modifiers
    logical, intent(out)   :: launch_done     !> `true` iff `rlaunch` read
    character(len=10), dimension(maxwords) :: dictionary !> keywords for parser
    integer,           dimension(maxwords) :: map        !> keycodes for parser
    integer, intent(out)   :: nwords          !> number of keywords
    integer, intent(in)    :: maxwords        !> dictionary & map capacity

    !! Local variables
    real, dimension(3,3) :: t
    real, dimension(3)   :: n1, n, nx, ny, n2, nz
    character(len=13)    :: down
    logical              :: nobuffer, atend
    character(len=80)    :: buffer, string, command
    real(dp_k), dimension(12,2) :: download

    integer :: i, j, k
    integer :: ntype, ntype2
    real    :: r1, r2

    open(unit=nbod, file=fbod, status='unknown')

    launch_done = .false.
    nelts = 0
    tolset = .false.
    unitcode = 'L '             ! unspecified length unit
    bt = .false.
    bm = .false.
    bw = .false.
    bc = .false.
    bbf = .false.
    hscale = 1.0             ! Set to unity, unless modified by HUNITS

    call webster(dictionary, map, nwords, maxwords)

    nobuffer = .true.
    call nextstring(nobuffer, buffer, string, atend)

    do while (.not. atend) !  GRANDDADDY loop: process command in file

       call testword(string, ntype, dictionary, map, nwords, maxwords)
       command = string

       if (ntype == sphere_code .or. ntype == cube_code) then

          down = 'rrrr@'
          call modifiers(nobuffer, buffer, string, atend, down, download, &
               &         dictionary, map, maxwords, nwords)
          nelts = nelts + 1
          eltype(nelts) = ntype
          do i = 1,4
             bv(nelts,i) = sngl(download(i,1))
          end do

       else if (ntype == triangle_code) then

          down = 'rrrrrrrrr@'
          call modifiers(nobuffer, buffer, string, atend, down, download, &
               &         dictionary, map, maxwords, nwords)
          nelts = nelts + 1
          eltype(nelts) = ntype
          do i = 1,9
             bv(nelts,i) = sngl(download(i,1))
          end do

       else if (ntype == donut_code) then

          down = 'rrrrrrrr@'
          call modifiers(nobuffer, buffer, string, atend, down, download, &
               &         dictionary, map, maxwords, nwords)
          nelts = nelts + 1
          eltype(nelts) = ntype
          do i = 1,8
             bv(nelts,i) = sngl(download(i,1))
          end do
          r1 = bv(nelts,7)
          r2 = bv(nelts,8)
          if (r2 > r1) then
             write (nzno,9009)
9009         format('Impossible torus:  r2 > r1 ')
             stop 'Impossible torus:  r2 > r1 '
          end if
          do j = 1,3
             n1(j) = bv(nelts,j+3)
          end do
          call normalize(n1,n)
          call zeerot(n,t)
          do j = 1,3
             do k = 1,3
                rotations(nelts,j,k) = t(j,k)
             end do
          end do
          do j = 1,3
             bv(nelts,j+3) = n(j)
          end do

       else if (ntype == disk_code) then

          down = 'rrrrrrr@'
          call modifiers(nobuffer, buffer, string, atend, down, download, &
               &         dictionary, map, maxwords, nwords)
          nelts = nelts + 1
          eltype(nelts) = ntype
          do i = 1,7
             bv(nelts,i) = sngl(download(i,1))
          end do
          do j = 1,3
             n1(j) = bv(nelts,j+3)
          end do
          call normalize(n1,n)
          call zeerot(n,t)
          do j = 1,3
             do k = 1,3
                rotations(nelts,j,k) = t(j,k)
             end do
          end do
          do j = 1,3
             bv(nelts,j+3) = n(j)
          end do

       else if (ntype == open_cylinder_code .or. &
            &   ntype == solid_cylinder_code) then

          down = 'rrrrrrrr@'
          call modifiers(nobuffer, buffer, string, atend, down, download, &
               &         dictionary, map, maxwords, nwords)
          nelts = nelts + 1
          eltype(nelts) = ntype
          do i = 1,8
             bv(nelts,i) = sngl(download(i,1))
          end do
          do j = 1,3
             n1(j) = bv(nelts,j+3)
          end do
          call normalize(n1,n)
          call zeerot(n,t)
          do j = 1,3
             do k = 1,3
                rotations(nelts,j,k) = t(j,k)
             end do
          end do
          do j = 1,3
             bv(nelts,j+3) = n(j)
          end do

       else if (ntype == ellipsoid_code) then
            
          down = 'rrrrrrrrrrrr@'
          call modifiers(nobuffer, buffer, string, atend, down, download, &
               &         dictionary, map, maxwords, nwords)
          nelts = nelts + 1
          eltype(nelts) = ntype
          do i = 1,12
             bv(nelts,i) = sngl(download(i,1))
          end do
          do j = 1,3
             n1(j) = bv(nelts,j+3)      
             n2(j) = bv(nelts,j+6)
          end do
          call normalize(n1,nx)
          call normalize(n2,ny)
          do j = 1,3
             bv(nelts,j+3) = nx(j)
             bv(nelts,j+6) = ny(j)
          end do
          call cross_product(nx,ny,nz)
          call xyzrot(nx,ny,nz,t)
          do j = 1,3
             do k = 1,3
                rotations(nelts,j,k) = t(j,k)
             end do
          end do

       else if (ntype == pillar_code) then

          down = 'rrrrrrrrrrrr@'
          call modifiers(nobuffer, buffer, string, atend, down, download, &
               &         dictionary, map, maxwords, nwords)
          nelts = nelts + 1
          eltype(nelts) = ntype
          do i = 1,12
             bv(nelts,i) = sngl(download(i,1))
          end do

       else if (ntype == skin_code) then
        
          down = 'r@'
          call modifiers(nobuffer,buffer, string, atend, down, download, &
               &         dictionary, map, maxwords, nwords)
          tol = sngl(download(1,1))
          tolset = .true.

       else if (ntype == bf_code) then

          down = 'd@'
          call modifiers(nobuffer, buffer, string, atend, down, download, &
               &         dictionary, map, maxwords, nwords)
          buoy(1) = download(1,1)
          buoy(2) = download(1,2)
          bbf = .true.

       else if (ntype == rlaunch_code) then
        
          down = 'r@'
          call modifiers(nobuffer, buffer, string, atend, down, download, &
               &         dictionary, map, maxwords, nwords)
          rlaunch = sngl(download(1,1))
          launch_done = .true.

       else if (ntype == units_code) then

          down = 'n@'
          call modifiers(nobuffer, buffer, string, atend, down, download, &
               &         dictionary, map, maxwords, nwords)
          ntype2 = nint(download(1,1))
          if (ntype2 == meter_code) then
             unitcode = 'm '
          else if (ntype2 == cm_code) then
             unitcode = 'cm'
          else if (ntype2 == nm_code) then
             unitcode = 'nm'
          else if (ntype2 == angstrom_code) then
             unitcode = 'A '
          else if (ntype2 == length_code) then
             unitcode = 'L '
          else 
             write(nzno,902)
902          format('Bad modifier to UNITS command')
             stop 'bad modifier to UNITS command'
          end if

       else if (ntype == hunits_code) then

          down = 'rn@'
          call modifiers(nobuffer, buffer, string, atend, down, download, &
               &         dictionary, map, maxwords, nwords)
          ntype2 = nint(download(2,1))
          hscale = sngl(download(1,1))
          if (ntype2 == meter_code) then
             unitcode = 'm '
          else if (ntype2 == cm_code) then
             unitcode = 'cm'
          else if (ntype2 == nm_code) then
             unitcode = 'nm'
          else if (ntype2 == angstrom_code) then
             unitcode = 'A '
          else 
             write(nzno,903)
903          format('Bad modifier to HUNITS command')
             stop 'bad modifier to HUNITS command'
          end if

       else if (ntype == temp_code) then

          down = 'dn@'
          call modifiers(nobuffer, buffer, string, atend, down, download, &
               &         dictionary, map, maxwords, nwords)
          temp(1) = download(1,1)
          temp(2) = download(1,2)
          ntype2 = nint(download(2,1))
          if (ntype2 == kelvin_code) then
             tunit = 'K     '
          else if (ntype2 == celcius_code) then
             tunit = 'C     '
          else 
             write(nzno,904)
904          format('Bad modifier to TEMP command')
             stop 'bad modifier to TEMP command'
          end if
          bt = .true.

       else if (ntype == mass_code) then

          down = 'dn@'
          call modifiers(nobuffer, buffer, string, atend, down, download, &
               &         dictionary, map, maxwords, nwords)
          mass(1) = download(1,1)
          mass(2) = download(1,2)
          ntype2 = nint(download(2,1))
          if (ntype2 == da_code) then
             munit = 'Da    '
          else if (ntype2 == kda_code) then
             munit = 'kDa   '
          else if (ntype2 == gram_code) then
             munit = 'g     '
          else if (ntype2 == kg_code) then
             munit = 'kg    '
          else 
             write(nzno,905)
905          format('Bad modifier to MASS command')
             stop 'bad modifier to MASS command'
          end if
          bm = .true.

       else if (ntype == visc_code) then
          down = 'dn@'
          call modifiers(nobuffer, buffer, string, atend, down, download, &
               &         dictionary, map, maxwords, nwords)
          if (.not. bw) then
             visc(1) = download(1,1)
             visc(2) = download(1,2)
             ntype2 = nint(download(2,1))
             if (ntype2 == poise_code) then
                vunit = 'p     '
             else if (ntype2 == cp_code) then
                vunit = 'cp    '
             else
                write(nzno,906)
906             format('Bad modifier to VISCOSITY command')
                stop 'bad modifier to VISCOSITY command'
             end if
             bc = .true.
          else
             write(nzno,711)
             write(nzno,712)
          end if
711       format('VISCOSITY and SOLVENT commands both found.')
712       format('Only the first was used.')

       else if (ntype == solvent_code) then
          down = 'n@'
          call modifiers(nobuffer, buffer, string, atend, down, download, &
               &         dictionary, map, maxwords, nwords)
          if (.not. bc) then
             ntype2 = nint(download(1,1))
             if (ntype2 /= water_code) then
                write(nzno,907)
907             format('Bad modifier to SOLVENT command')
                stop 'bad modifier to SOLVENT command'
             end if
              bw = .true.
          else
             write(nzno,711)
             write(nzno,712)
          end if

       else
          write(nzno,610) 
          write(nzno,611) command
610       format('Unrecognized command: ')
611       format(a80)
          stop 'unrecognized command'

       end if

       if (nelts > maxelts) then
          write(nzno,612)
612       format('Too many elements.')
          stop 'too many elements'
       end if

       call nextstring(nobuffer, buffer, string, atend)

    end do              !  end of GRANDDADDY loop

    if (nelts == 0) then
       write(nzno,613)
613    format('No elements loaded')
       stop 'no elements loaded'
    end if

    close(nbod)

  end subroutine parse

  !! ================================================================
  !!
  !! Private routines
  !!
  !! ================================================================

  !> @brief builds an array of recognized keywords and the corresponding
  !> array of integer keycodes.
  !>
  !> @warning should be much simpler than this.
  !> 
  !> @todo replace by data structure that is \e statically initialized
  !> to be alphabetically sorted so it can be searched using a binary
  !> search routine.

  subroutine webster(dictionary, map, nwords, maxwords)

    use zeno_codes_data

    !! argument declarations

    character(len=10), dimension(maxwords) :: dictionary !> keywords list
    integer,           dimension(maxwords) :: map        !> keycodes list
    integer, intent(out) :: nwords                       !> # of keywords
    integer, intent(in)  :: maxwords                     !> list capacity

    !! local variables

    character(len=10) :: word
    character(len=51) :: &
         & slash1, slash2, slash3, slash4, slash5, slash6, slash7, &
         & copy

    integer :: i

    nwords = 65

    slash1 = 'ELLIPSOID,ellipsoid,VISCOSITY,viscosity,TRIANGLE,  '
    copy = slash1
    do i = 1, 5
       call commaparse(copy, word)
       dictionary(i) = word
    end do
    map(1) = ellipsoid_code
    map(2) = ellipsoid_code
    map(3) = visc_code
    map(4) = visc_code
    map(5) = triangle_code

    slash2 = 'triangle,CYLINDER,cylinder,SOLVENT,solvent,RLAUNCH,'
    copy = slash2
    do i = 6, 11
       call commaparse(copy, word)
       dictionary(i) = word
    end do
    map(6) = triangle_code
    map(7) = open_cylinder_code
    map(8) = open_cylinder_code
    map(9) = solvent_code
    map(10) = solvent_code
    map(11) = rlaunch_code

    slash3 = 'rlaunch,SOL_CYL,sol_cyl,SPHERE,sphere,HUNITS,      '
    copy = slash3
    do i = 12, 17
       call commaparse(copy, word)
       dictionary(i) = word
    end do
    map(12) = rlaunch_code
    map(13) = solid_cylinder_code
    map(14) = solid_cylinder_code
    map(15) = sphere_code
    map(16) = sphere_code
    map(17) = hunits_code

    slash4 = 'hunits,PILLAR,pillar,UNITS,units,TORUS,torus,WATER,' 
    copy = slash4
    do i = 18, 25
       call commaparse(copy, word)
       dictionary(i) = word
    end do
    map(18) = hunits_code
    map(19) = pillar_code
    map(20) = pillar_code
    map(21) = units_code
    map(22) = units_code
    map(23) = donut_code
    map(24) = donut_code
    map(25) = water_code

    slash5 = 'water,DISK,disk,CUBE,cube,MASS,mass,TEMP,temp,kDa, '
    copy = slash5
    do i = 26, 35
       call commaparse(copy,word)
       dictionary(i) = word
    end do
    map(26) = water_code
    map(27) = disk_code
    map(28) = disk_code
    map(29) = cube_code
    map(30) = cube_code
    map(31) = mass_code
    map(32) = mass_code
    map(33) = temp_code
    map(34) = temp_code
    map(35) = kda_code

    slash6 = 'TO,to,ST,st,SC,sc,cm,nm,Da,kg,cp,PI,pi,S,s,T,t,D,d,'
    copy = slash6
    do i = 36, 54
       call commaparse(copy,word)
       dictionary(i) = word
    end do
    map(36) = donut_code
    map(37) = donut_code
    map(38) = skin_code
    map(39) = skin_code
    map(40) = solid_cylinder_code
    map(41) = solid_cylinder_code
    map(42) = cm_code
    map(43) = nm_code
    map(44) = da_code
    map(45) = kg_code
    map(46) = cp_code
    map(47) = pillar_code
    map(48) = pillar_code
    map(49) = sphere_code
    map(50) = sphere_code
    map(51) = triangle_code
    map(52) = triangle_code
    map(53) = disk_code
    map(54) = disk_code

    slash7 = 'C,E,e,m,A,L,K,g,p,BF,bf,                           '
    copy = slash7
    do i = 55, 65
       call commaparse(copy, word)
       dictionary(i) = word
    end do
    map(55) = celcius_code
    map(56) = ellipsoid_code
    map(57) = ellipsoid_code
    map(58) = meter_code
    map(59) = angstrom_code
    map(60) = length_code
    map(61) = kelvin_code
    map(62) = gram_code
    map(63) = poise_code
    map(64) = bf_code
    map(65) = bf_code

  end subroutine webster

  !! ----------------------------------------------------------------

  !> Uses commas to tokenize a string
  !!
  !! @warning modifies `copy` argument

  subroutine commaparse(copy, word)

    !! argument declarations

    character(len=51) :: copy
    character(len=10) :: word

    !! local variables

    character(len=51) :: list
    integer           :: ii, nzip

    nzip = index(copy, ',')
    word(1:nzip-1) = copy(1:nzip-1)
    word(nzip:nzip) = '@'
    if (nzip < 10) then
       do ii = nzip+1,10
          word(ii:ii) = ' '
       end do
    end if

    list(1:51-nzip) = copy(nzip+1:51)
    do ii = 52-nzip,51
       list(ii:ii) = ' '
    end do
    copy = list

  end subroutine commaparse

  !! ----------------------------------------------------------------

  !> @brief Tests if word is a recognized keyword
  !!
  !! Is string a recognized word?  If so encode its value in `ntype`
  !! Here is the complete dictionary: if string does not match one of
  !! these, abort, otherwise return its equivalent code.
  !! 
  !! @verbatim
  !!       Word              Code
  !!----------------------------------
  !!    ELLIPSOID       ellipsoid_code
  !!    ellipsoid       ellipsoid_code
  !!    VISCOSITY       visc_code
  !!    viscosity       visc_code
  !!    TRIANGLE        triangle_code
  !!    triangle        triangle_code
  !!    CYLINDER        open_cylinder_code
  !!    cylinder        open_cylinder_code
  !!    SOLVENT         solvent_code
  !!    solvent         solvent_code
  !!    RLAUNCH         rlaunch_code
  !!    rlaunch         rlaunch_code
  !!    SOL_CYL         solid_cylinder_code
  !!    sol_cyl         solid_cylinder_code
  !!    SPHERE          sphere_code
  !!    sphere          sphere_code
  !!    HUNITS          hunits_code
  !!    hunits          hunits_code
  !!    PILLAR          pillar_code
  !!    pillar          pillar_code
  !!    UNITS           units_code
  !!    units           units_code
  !!    TORUS           donut_code
  !!    torus           donut_code
  !!    WATER           water_code
  !!    water           water_code
  !!    DISK            disk_code
  !!    disk            disk_code
  !!    CUBE            cube_code
  !!    cube            cube_code
  !!    MASS            mass_code
  !!    mass            mass_code
  !!    TEMP            temp_code
  !!    temp            temp_code
  !!    kDa             kda_code
  !!    TO              donut_code
  !!    to              donut_code
  !!    ST              skin_code
  !!    st              skin_code
  !!    SC              solid_cylinder_code
  !!    sc              solid_cylinder_code
  !!    cm              cm_code
  !!    nm              nm_code
  !!    Da              da_code
  !!    kg              kg_code
  !!    cp              cp_code
  !!    PI              pillar_code
  !!    pi              pillar_code
  !!    S               sphere_code
  !!    s               sphere_code
  !!    T               triangle_code
  !!    t               triangle_code
  !!    D               disk_code
  !!    d               disk_code
  !!    C               celcius_code
  !!    E               ellipsoid_code
  !!    e               ellipsoid_code
  !!    m               meter_code
  !!    A               angstrom_code
  !!    L               length_code
  !!    K               kelvin_code
  !!    g               gram_code
  !!    p               poise_code
  !!    BF              bf_code
  !!    bf              bf_code
  !! @endverbatim
  !!
  !! @todo
  !! Replace the whole thing by a more suitable data structure (e.g., a
  !! hash table or simple enums generated by a lexical analyzer)

  subroutine testword(string, ntype, dictionary, map, nwords,&
       & maxwords)

    use zeno_codes_data
    use zeno_files_data

    !! declare arguments

    character(len=80) :: string
    integer           :: ntype
    character(len=10), dimension(maxwords) :: dictionary !> keywords list
    integer,           dimension(maxwords) :: map        !> keycodes list
    integer, intent(in)                    :: nwords     !> # of keywords
    integer, intent(in)                    :: maxwords   !> list capacity

    !! local variables

    integer :: i, nlen, match

    nlen = index(string,' ') - 1

    if (nlen > 9) then
       write(nzno, 900)
       write(nzno,901) string
900    format('Word longer than any in dictionary:  ')
901    format(a80)
       stop 'word longer than any in dictionary'
    end if

    match = 0
    do i = 1,nwords
       if (dictionary(i)(nlen+1:nlen+1) == '@') then
          if (dictionary(i)(1:nlen) == string(1:nlen)) then
             match = i
          end if
       end if
    end do

    if (match == 0) then
       write(nzno, 800)
       write(nzno,901) string
800    format('Word not in dictionary:  ')
       stop 'word not in dictionary'
    end if

    ntype = map(match)

  end subroutine testword

  !! ----------------------------------------------------------------

  !> Interprets the next input string

  subroutine nextstring(nobuffer, buffer, string, atend)

    use zeno_files_data

    !! declare arguments

    logical              :: nobuffer
    character(len=80)    :: buffer
    character(len=80)    :: string
    logical, intent(out) :: atend

    !! local variables

    integer :: i, next

    atend = .false.

    do
       if (nobuffer) then
          read(nbod,980,end=20) buffer
          call cleartabs(buffer)

          !> if (.not. silent) print 980,buffer
          nobuffer = .false.
980       format(a80)
       end if

       if (buffer(1:1) == '*') then
          nobuffer = .true.
          cycle
       end if

       call packspaces(buffer,next)
       if (next == 0) then
          nobuffer = .true.
          cycle
       end if

       exit
    end do

    do i = 1,80
       string(i:i) = ' '
    end do
    string(1:next-1) = buffer(1:next-1)
    do i = 1,next-1
       buffer(i:i) = ' '
    end do

    return

20  continue
    atend = .true.
    return

  end subroutine nextstring

  !! ----------------------------------------------------------------

  !> Replaces `tab` characters by `blankspaces`
 
  subroutine cleartabs(buffer)

    !! declare arguments

    character(len=80) :: buffer

    integer :: i

    !! Convert tab-characters in input buffer to spaces

    do i = 1,80
       if (ichar(buffer(i:i)) == 9) buffer(i:i) = ' '
    end do
  
  end subroutine cleartabs

  !! ----------------------------------------------------------------

  !> Pack leading spaces from the string `buffer`
  !! 
  !! @warning Does too much copying.  Line parser should simply discard
  !! leading white space!
 
  subroutine packspaces(buffer, next)

    !! declare arguments

    character(len=80) :: buffer, copy

    !! Local variables

    integer :: next, ngo

    ngo = 0

    do while (buffer(1:1) == ' ')
       copy = buffer
       buffer(1:79) = copy(2:80)
       buffer(80:80) = ' '
       ngo = ngo + 1
       if (ngo == 80) then
          next = 0
          return
       end if
    end do

    next = index(buffer, ' ')

  end subroutine packspaces

  !! ----------------------------------------------------------------

  !> Modifies commad read at beginning of line
  !! 
  !! @todo Replace by a grammar driven approach

  subroutine modifiers(nobuffer, buffer, string, atend, down, &
       &     download, dictionary, map, maxwords, nwords)

    use zeno_files_data

    !! Declare arguments

    logical :: nobuffer         !> empty buffer flag
    character(len=80) :: buffer !> rest of line
    character(len=80) :: string !> copy of rest of line?
    logical :: atend            !> reached end of line?
    character(len=13) :: down   !> unknwon?
    real(dp_k), dimension(12,2) :: download              !> numbers ??
    character(len=10), dimension(maxwords) :: dictionary !> known keywords
    integer,           dimension(maxwords) :: map        !> keycodes list
    integer, intent(in)                    :: nwords     !> # keywords
    integer, intent(in)                    :: maxwords   !> list capacity

    !! Local variables

    real(dp_k), dimension(2) :: w
    integer :: i, ntype

    do i = 1,13

       if (down(i:i) == '@') then

          return

       else if (down(i:i) == 'n') then

          call nextstring(nobuffer, buffer, string, atend)

          if (atend) then
             write(nzno,400)
400          format('Unexpected end of body file')
             stop 'unexpected end of body file'
          end if
              
          call testword(string, ntype, dictionary, map, nwords, maxwords)
          download(i,1) = dble(ntype)
          download(i,2) = 0.0d0

       else if (down(i:i) == 'r') then

          call nextstring(nobuffer, buffer, string, atend)
            
          if (atend) then
             write(nzno,400)
             stop 'unexpected end of body file'
          end if

          read(string,*,err=999) download(i,1)
          download(i,2) = 0.0d0

       else if (down(i:i) == 'd') then

          call nextstring(nobuffer,buffer,string,atend)
            
          if (atend) then
             write(nzno,400)
             stop 'unexpected end of body file'
          end if

          call floatstring(string,w)

          download(i,1) = w(1)
          download(i,2) = w(2)

       end if
    end do
    return

999 continue

    write(nzno, 620)
620 format('Trying to interpret this string as a number: ')
    write(nzno,621) string
621 format(a80)
    stop 'number conversion error'

  end subroutine modifiers

  !! ----------------------------------------------------------------

  !> Removes leading white space
  !! 
  !! @warning Does a lot of copying!
  !! 
  !! @todo Either replace or count leading spaces and then copy

  subroutine pack80(line)

    !! Declare arguments

    character(len=80) :: line

    !! Local variables

    character(len=80) :: copy
    integer           :: nz

    nz = 0

    do while (line(1:1) == ' ' .and. nz < 90)
       copy(1:79) = line(2:80)
       copy(80:80) = ' '
       line = copy
       nz = nz + 1
    end do

  end subroutine pack80

  !! ----------------------------------------------------------------

  !> parses and validates a `float` in `tline`
  !! 
  !! @return double value in `t1`
  !! 
  !! @warning should take advantage of facilities in runtime library!

  subroutine floatstring(tline, t1)

    !! Declare arguments

    character(len=80) :: tline
    real(dp_k), dimension(2) :: t1

    !! Local variables

    character(len=80) :: cline
    integer  :: i, j1, j2, j3
    integer  :: ndot, nex, nj
    real(dp_k) :: pre, eval, val, dig, place

    j1 = index(tline,'E')
    if (j1 == 0) then
       j2 = index(tline,'e')
       if (j2 == 0) then
          j3 = index(tline,' ')
          tline(j3:j3+1) = 'e0'
       end if
    else
       tline(j1:j1) = 'e'
    end if

    if (tline(1:1) == '-') then
       pre = -1.0d0
       tline(1:1) = ' '
    else if (tline(1:1) == '+') then
       pre = +1.0d0
       tline(1:1) = ' '
    else
       pre = 1.0d0
    end if

    call pack80(tline)

    j1 = index(tline,'e')
    if (j1 == 0) stop 'e should be there by now'

    do i = 1,j1
       cline(i:i) = ' '
    end do
    do i = j1+1,80
       cline(i:i) = tline(i:i)
       tline(i:i) = ' '
    end do
    tline(j1:j1)= ' '

    read(cline,*) nex

    nj = index(tline,'.')
    if (nj == 0) then
       nj = index(tline,' ')
       tline(nj:nj) = '.'
    end if

    ndot = index(tline,'.')
    if (ndot == 0) stop 'dot should be there by now'

    nj = index(tline,' ') - 1
        
    if (nj > ndot) then
       eval = 10.0d0**(ndot-nj)
    else if (nj == ndot) then
       eval = 1.0d0
    else
       stop 'dotty'
    end if

    val = 0.0d0
    do i = 1,nj
       if (i /= ndot) then
          read(tline(i:i),*) dig
          if (i < ndot) then
             place = 10.0d0 ** (ndot-i-1)
          else 
             place = 10.0d0 ** (ndot-i)
          end if
          val = val + dig*place
       end if
    end do

    t1(1) = pre*val*(10.0d0**nex)
    t1(2) = 0.5d0*eval*(10.0d0**nex)

  end subroutine floatstring

  !! ================================================================

end module zeno_parser

!!! ================================================================

!!! Local Variables:
!!! mode: f90
!!! time-stamp-line-limit: 30
!!! fill-column: 72
!!! End:
