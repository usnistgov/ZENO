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
!! This file provides the `zeno_setup` module which defines routines
!! invoked as part of the initialization of ZENO.  Among other things,
!! these routines handle command line options.
!! 
!! @author Mark Mansfield (formerly of Stevens Institute) and students
!! @author Walid Keyrouz (walid.keyrouz@nist.gov)
!! @date   Tue Oct 22 15:59:15 2013 EDT
!!
!! @defgroup zeno_setup Initialization Routines
!! 
!! Groups together the subroutines that are used in setting up ZENO.
!! These routines are typically called when the program first executes.
!
! Time-stamp: <2015-01-05 16:44:50 walid>
!
! ================================================================

module zeno_setup

  !! ================================================================

  use getopt_m
  use zeno_enums
  use zeno_debug,   only : silent, print_flush
  use zeno_options, only : eig_opt, rng_opt
  use zeno_rng,     only : seeder
  use zeno_output,  only : gettime
  use zeno_version, only : VERSION_STRING

  implicit none
  private

  public setup

contains

  !! ================================================================
  !!
  !! Public routines
  !!
  !! ================================================================

  !> Parse the invocation string, and initialize random numbers
  !! 
  !! @todo rename to reflect functionality
  !!
  !! @warning argument list is too long.

  subroutine setup(id, m1, actions, start, actout, &
       &           savehits_z, savehits_i, savehits_s, savehits_p)

    use zeno_files_data

    !! Declare arguments
    character(len=25)     :: id        !> ZENO case id
    integer, dimension(4) :: m1        !> Num of MC steps in computation i
    character(len=4)      :: actions   !> TBD
    character(len=28)     :: start     !> TBD
    character(len=30), dimension(4) ::  actout !> TBD
    logical :: savehits_z              !> TBD
    logical :: savehits_i              !> TBD
    logical :: savehits_s              !> TBD
    logical :: savehits_p              !> TBD

    !! Local variables
    integer :: i
    integer :: mult
    integer :: nac, nfo, nsp_id, ntg, nzip
    character(len=30) :: ac
    character(len=30) :: pix
    character(len=32) :: com
    character(len=25) :: dfl

    !! Options processing
    type(option_s):: zeno_opts(8)

    character(len=25) :: cli_id    !> CLI case id
    character(len=10) :: cli_od    !> CLI output directory
    character(len=28) :: cli_ts    !> CLI initial seed
    logical :: flag_od, flag_ts    !> flags indicating CLI presence
    logical :: flag_id
    logical :: flag_version
    integer :: len_od, nsp_od      !> Output dir strlen & non-blank idx
    integer :: len_id, nsp_oid
    character(len=10) :: cli_eig   !> CLI eigen choice
    integer :: len_eig             !> eigen choice length
    character(len=10) :: cli_rng   !> CLI rng choice
    integer :: len_rng             !> RNG choice length

    character(len=25) :: dir_id    !> Augmented pathname: odir/id

    flag_id = .false.
    flag_od = .false.
    flag_ts = .false.

    flag_version = .false.

    savehits_z = .false.
    savehits_i = .false.
    savehits_s = .false.
    savehits_p = .false.

    !> Process the following command-line options:
    !! 
    !! - `--case_id` or `-c` -- Case ID
    !! - `--seed` or `-s` -- Seed of RNG
    !! - `--odir` or `-o` -- Output directory
    !! - `--version` or `-V` -- Version number
    !! - `--verbose` or `-v` -- Verbose output
    !! - `--flush_print` or `-f` -- Flush print to show progress output
    !!   (similar to `curses`)
    !! - `--eigen` or `-e` -- Choose Eigenvalue implementation
    !! - `--rng` or `-r` -- Choose RNG implementation
    !! 
    !! using the `getopt_m:getopt` function.
    !! 
    !! These options must appear at the beginning of the command line
    !! because this version of getopt does not rearrange the command
    !! line arguments as is the case in its C/C++ counterpart.
    !!
    !! getopt.f90 can be retrieved from:
    !! http://lagrange.mechse.illinois.edu/mwest/partmc/partmc-2.2.0/src
    !!
    !! Please see the copyright at the beginning of the file.

    zeno_opts(1) = option_s( "seed",        .true.,  's' )
    zeno_opts(2) = option_s( "odir",        .true.,  'd' )
    zeno_opts(3) = option_s( "case_id",     .true.,  'c' )
    zeno_opts(4) = option_s( "version",     .false., 'V' )
    zeno_opts(5) = option_s( "verbose",     .false., 'v' )
    zeno_opts(6) = option_s( "flush_print", .false., 'f' )
    zeno_opts(7) = option_s( "eigen",       .true.,  'e' )
    zeno_opts(8) = option_s( "rng",         .true.,  'r' )

    cli_opts_lbl: do
       select case( getopt( "s:d:c:vfe:r:", zeno_opts ))

       case( char(0))
          exit

       case ( 'c' )             !> `case_id`: basename of output files.
          flag_id   = .true.
          len_id    = len_trim(optarg)
          cli_id(:) = optarg(:len_id)   !! strlen(case_id) <= 25

       case( 'd' )              !> `odir`: output dir name
          flag_od   = .true.
          len_od    = len_trim(optarg)
          cli_od(:) = optarg(:len_od)   !! strlen(odir) <= 6

       case( 's' )              !> `seed`: date string used as seed
          flag_ts   = .true.
          cli_ts(:) = optarg(:28)

       case( 'V' )              !> `version`: output version number
          flag_version = .true.
          write(*,"((a)(a))") 'ZENO version: ', VERSION_STRING
          stop

       case( 'v' )              !> `verbose`: debug output
          silent = .false.

       case( 'f' )              !> `flush_print`: shows progress via BS char
          print_flush = .true.

       case ( 'e' )             !> `eigen`: eigenvalue routines used.
          !> Only Wikipedia for now
          len_eig    = len_trim(optarg)
          cli_eig(:) = optarg(:len_eig)

          if (cli_eig(1:1) == 'W' .or. cli_eig(1:1) == 'w') then
             eig_opt = EIG_wp   !! Use eigenvalue routine from Wikipedia
          else
             eig_opt = EIG_wp
             write(*,*) 'Warning: illegal eigen method specified--', cli_eig
          end if

       case( 'r' )              !> `RNG` routine library used.
          !> Only GSL (FGSL) for now
          len_rng = len_trim(optarg)
          cli_rng(:) = optarg(:len_rng)

          if (cli_rng(1:1) == 'g' .or. cli_rng(1:1) == 'G') then
             rng_opt = RNG_fgsl !! PRNG from [f]GSL
          else
             rng_opt = RNG_fgsl
             write(*,*) 'Warning: illegal RNG method specified--', cli_rng
          end if

       case( '?' )
          print '(a,a)', 'unknown option: ', optopt
          stop

       case default
          print '(a,a,a)', 'unhandled option: ', optopt, ' (this is a bug)'
       end select

    end do cli_opts_lbl

    !! Get the body ID:
    !! Index optind names body file.
    call getarg(optind,id)

    !! Get the action codes
    ntg = iargc()
    nac = 0

    !! Parse command line after ".bod" file.
    do i = optind+1,ntg
       nfo = 0
       call getarg(i,ac)

       if (ac(1:1) == 'z') then

          nac = nac + 1
          actions(nac:nac) = 'z'
          nfo = 1

       else if (ac(1:1) == 'Z') then

          nac = nac + 1
          actions(nac:nac) = 'Z'
          savehits_z = .true.
          nfo = 1

       else if (ac(1:1) == 'i') then

          nac = nac + 1
          actions(nac:nac) = 'i'
          nfo = 1

       else if (ac(1:1) == 'I') then

          nac = nac + 1
          actions(nac:nac) = 'I'
          savehits_i = .true.
          nfo = 1

       else if (ac(1:1) == 's') then

          nac = nac + 1
          actions(nac:nac) = 's'
          nfo = 1

       else if (ac(1:1) == 'S') then

          nac = nac + 1
          actions(nac:nac) = 'S'
          savehits_s = .true.
          nfo = 1

       else if (ac(1:1) == 'p') then

          nac = nac + 1
          actions(nac:nac) = 'p'
          nfo = 1

       else if (ac(1:1) == 'P') then

          nac = nac + 1
          actions(nac:nac) = 'P'
          nfo = 1
          savehits_p = .true.

       else if (ac(1:1) == 'c') then

          nac = nac + 1
          actions(nac:nac) = 'c'
          nfo = 1

       else if (ac(1:1) == 'C') then

          nac = nac + 1
          actions(nac:nac) = 'C'
          nfo = 1
          savehits_i = .true.

       else
          stop 'bad action code -- expecting:  zisp or ZISP'
       end if

       if (nac > 4) stop 'too many actions!!!!'

       if (nfo == 1) then

          actout(nac) = ac
          ac(1:1) = ' '
          pix = ac

          mult = 1
          nzip = index(pix,'t')
          if (nzip /= 0) then
             mult = 1000*mult
             pix(nzip:nzip) = ' '
          end if

          nzip = index(pix,'m')
          if (nzip /= 0) then
             mult = 1000000*mult
             pix(nzip:nzip) = ' '
          end if

          nzip = index(pix,'b')
          if (nzip /= 0) then
             mult = 1000000000*mult
             pix(nzip:nzip) = ' '
          end if

          read(pix,*) m1(nac)
          m1(nac) = m1(nac) * mult

       end if

    end do

    if (nac < 4) then
       do i = nac+1,4
          actions(i:i) = '.'
       end do
    end if

    if (.not. flag_id) then
       cli_id(:) = char(0)
       cli_id(1:len_trim(id)) = id(:)
    end if

    nsp_id  = index(id,' ')
    nsp_oid = index(cli_id, ' ')

    dir_id(:) = char(0)

    if (flag_od) then
       !! Warning: must fit in 25 characters!
       dir_id(1:len_od)          = cli_od(1:len_od)
       dir_id(len_od+1:len_od+1) = '/'
       nsp_od = len_od + 2
    else
       nsp_od = 1
    end if

    if (flag_id) then
       dir_id(nsp_od:nsp_od+nsp_oid-1) = cli_id(1:nsp_oid-1)
       nsp_od = nsp_od + nsp_oid - 1
    else
       dir_id(nsp_od:nsp_od+nsp_id-1)  = id(1:nsp_id-1)
       nsp_od = nsp_od + nsp_id - 1
    end if

    !! Assume 'bod' file in current dir.
    fbod = id
    fbod(nsp_id:nsp_id+3) = '.bod'

    !! Other files go into output directory
    fzno = dir_id
    fzno(nsp_od:nsp_od+3) = '.zno'
    fznr = dir_id
    fznr(nsp_od:nsp_od+3) = '.znr'
    fstk = dir_id
    fstk(nsp_od:nsp_od+3) = '.stk'
    fdfl = dir_id
    fdfl(nsp_od:nsp_od+3) = '.dfl'
    fefl = dir_id
    fefl(nsp_od:nsp_od+3) = '.efl'
    fzsd = dir_id
    fzsd(nsp_od:nsp_od+3) = '.zsd'
    fzh  = dir_id
    fzh (nsp_od:nsp_od+2) = '.zh'
    fih  = dir_id
    fih (nsp_od:nsp_od+2) = '.ih'
    fsh  = dir_id
    fsh (nsp_od:nsp_od+2) = '.sh'
    fph  = dir_id
    fph (nsp_od:nsp_od+2) = '.ph'
    flog = dir_id
    flog(nsp_od:nsp_od+3) = '.log'

    nbod = 20
    nzno = 99
    nznr = 98
    nstk = 91
    ndfl = 82
    nefl = 82
    nzsd = 83
    nzh = 25
    nih = 35
    nsh = 45
    nph = 55
    nlog = 84

    open(unit=nzno,file=fzno,status='unknown')
    open(unit=nlog,file=flog,status='unknown')

    com(1:7) = 'date > '
    !!          1234567
    com(8:32) = fdfl
    call system(com)

    dfl = fdfl
    call gettime(dfl,start)

    !! WK: handle "--seed" CLI option
    if (flag_ts) then
       call seeder(cli_ts)

       open(unit=nzsd,file=fzsd,status='unknown')
       rewind nzsd
       write(nzsd,'(A28)') cli_ts
       close(nzsd)
    else
       call seeder(start)
    end if

  end subroutine setup

  !! ================================================================

end module zeno_setup

!!! ================================================================

!!! Local Variables:
!!! mode: f90
!!! time-stamp-line-limit: 30
!!! fill-column: 72
!!! End:
