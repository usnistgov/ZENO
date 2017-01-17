! ------------------------------------------------------------
!
!> @file
!! Code released by Mark Gates, under GPL2 or later.  It provides the
!! getopt_m module which implements Fortran 95 functions, `getopt()`
!! and `getopt_long()`, with functionalities similar to those found in
!! GNU Libc.
!! @copyright 2008 by Mark Gates
!!
!! @defgroup getopt_m GetOpt Module
!! 
!! Released by Mark Gates under GPL2 or later, provides Fortran 95
!! functions, `getopt()` and `getopt_long()`, with functionalities
!! similar to those available in GNU Libc.
!! 
!! @details
!! @verbatim
!! ------------------------------------------------------------
!! Copyright 2008 by Mark Gates
!!
!! This program is free software; you can redistribute or modify it under
!! the terms of the GNU general public license (GPL), version 2 or later.
!!
!! This program is distributed in the hope that it will be useful, but
!! WITHOUT ANY WARRANTY; without even the implied warranty of
!! merchantability or fitness for a particular purpose.
!!
!! If you wish to incorporate this into non-GPL software, please contact
!! me (Mark Gates) regarding licensing terms.
!!
!! ------------------------------------------------------------
!! @endverbatim
!! The getopt_m module, released by Mark Gates under GPL2 or later,
!! provides Fortran 95 functions, `getopt()` and `getopt_long()`, with
!! functionalities similar to those available in GNU Libc.
!! - A copy of this file was retrieved by Walid Keyrouz
!!   (walid.keyrouz@nist.gov) from:
!!   http://lagrange.mechse.illinois.edu/mwest/partmc/partmc-2.2.0/src/getopt.F90
!! - It was edited only to enable documentation generation by Doxygen.
!! 
!! Fortran 95 `getopt()` and `getopt_long()`, similar to those in
!! standard C library.
!!
!!     ch = getopt( optstring, [longopts] )
!! 
!! Returns next option character from command line arguments.
!! If an option is not recognized, it returns '?'.
!! If no options are left, it returns a null character, char(0).
!! 
!! - `optstring` contains characters that are recognized as options.
!!   If a character is followed by a colon, then it takes a required
!!   argument.  For example, `"x"` recognizes `"-x"`, while `"x:"`
!!   recognizes `"-x arg"` or `"-xarg"`.
!!
!! - `optopt` is set to the option character, even if it isn't recognized.
!! 
!! - `optarg` is set to the option's argument.
!! 
!! - `optind` has the index of the next argument to process.  Initially
!!   `optind=1`.
!! 
!! Errors are printed by default. Set `opterr=.false.` to suppress them.
!! 
!! Grouped options are allowed, so @c "-abc" is the same as @c "-a -b -c".
!!
!! If `longopts` is present, it is an array of type(`option_s`), where
!! each entry describes one long option.
!!
!! ```
!!    type option_s
!!        character(len=80) :: name
!!        logical           :: has_arg
!!        character         :: val
!!    end type
!! ```
!!
!! The name field is the option name, without the leading @c -- double
!! dash.  Set the `has_arg` field to true if it requires an argument,
!! false if not.  The `val` field is returned.  Typically this is set
!! to the corresponding short option, so short and long options can be
!! processed together.  (But there is no requirement that every long
!! option has a short option, or vice-versa.)
!!
!! -----
!! 
!! EXAMPLE
!! @verbatim
!! program test
!!     use getopt_m
!!     implicit none
!!     character:: ch
!!     type(option_s):: opts(2)
!!     opts(1) = option_s( "alpha", .false., 'a' )
!!     opts(2) = option_s( "beta",  .true.,  'b' )
!!     do
!!         select case( getopt( "ab:c", opts ))
!!             case( char(0))
!!                 exit
!!             case( 'a' )
!!                 print *, 'option alpha/a'
!!             case( 'b' )
!!                 print *, 'option beta/b=', optarg
!!             case( '?' )
!!                 print *, 'unknown option ', optopt
!!                 stop
!!             case default
!!                 print *, 'unhandled option ', optopt, ' (this is a bug)'
!!         end select
!!     end do
!! end program test
!! @endverbatim
!!
!! Differences from C version:
!! - when options are finished, C version returns @c -1 instead of @c
!!   char(0), and thus stupidly requires an `int` instead of a `char`.
!! - does not support `optreset`
!! - does not support `"--"` as last argument
!! - if no argument, `optarg` is blank, not `NULL`
!! - `argc` and `argv` are implicit
!!
!! Differences for long options:
!! - optional argument to `getopt()`, rather than separate function
!!   `getopt_long()`
!! - `has_arg` is logical, and does not support `optional_argument`
!! - does not support flag field (and thus always returns val)
!! - does not support `longindex`
!! - does not support `"--opt=value"` syntax, only `"--opt value"`
!! - knows the length of `longopts`, so does not need an empty last record
!! 
!! @date   Last access of original file: 2014-12-22
!! @author Mark Gates
!! @author Modified by Walid Keyrouz (walid.keyrouz@nist.gov)
!!

module getopt_m
  implicit none

  character(len=80):: optarg
  character:: optopt
  integer:: optind=1
  logical:: opterr=.true.

  type option_s
     character(len=80) :: name
     logical           :: has_arg
     character         :: val
  end type option_s

  ! grpind is index of next option within group, always >= 2
  integer, private:: grpind=2

  public getopt

  private substr, process_long, process_short

contains

  ! ----------------------------------------
  ! Return str(i:j) if 1 <= i <= j <= len(str),
  ! else return empty string.
  ! This is needed because Fortran standard allows but doesn't *require* short-circuited
  ! logical AND and OR operators. So this sometimes fails:
  !     if ( i < len(str) .and. str(i+1:i+1) == ':' ) then
  ! but this works:
  !     if ( substr(str, i+1, i+1) == ':' ) then

  character function substr( str, i, j )
    ! arguments
    character(len=*), intent(in):: str
    integer, intent(in):: i, j

    if ( 1 <= i .and. i <= j .and. j <= len(str)) then
       substr = str(i:j)
    else
       substr = ''
    endif
  end function substr


  ! ----------------------------------------

  !> Function `getopt`
  !! @return `character`
  !! @param optstring input string specifying short options
  !! @param longopts  *optional* input string specifying long options
  !<
  character function getopt( optstring, longopts )
    ! arguments
    character(len=*), intent(in):: optstring
    type(option_s),   intent(in), optional:: longopts(:)

    ! local variables
    character(len=80):: arg

    optarg = ''
    if ( optind > iargc()) then
       getopt = char(0)
    endif

    call getarg( optind, arg )
    if ( present( longopts ) .and. arg(1:2) == '--' ) then
       getopt = process_long( longopts, arg )
    elseif ( arg(1:1) == '-' ) then
       getopt = process_short( optstring, arg )
    else
       getopt = char(0)
    endif
  end function getopt


  ! ----------------------------------------
  character function process_long( longopts, arg )
    ! arguments
    type(option_s),   intent(in):: longopts(:)
    character(len=*), intent(in):: arg

    ! local variables
    integer:: i

    ! search for matching long option
    optind = optind + 1
    do i = 1, size(longopts)
       if ( arg(3:) == longopts(i)%name ) then
          optopt = longopts(i)%val
          process_long = optopt
          if ( longopts(i)%has_arg ) then
             if ( optind <= iargc()) then
                call getarg( optind, optarg )
                optind = optind + 1
             elseif ( opterr ) then
                print '(a,a,a)', "Error: option '", trim(arg), "' requires an argument"
             endif
          endif
          return
       endif
    end do
    ! else not found
    process_long = '?'
    if ( opterr ) then
       print '(a,a,a)', "Error: unrecognized option '", trim(arg), "'"
    endif
  end function process_long


  ! ----------------------------------------
  character function process_short( optstring, arg )
    ! arguments
    character(len=*), intent(in):: optstring, arg

    ! local variables
    integer:: i, arglen

    arglen = len( trim( arg ))
    optopt = arg(grpind:grpind)
    process_short = optopt

    i = index( optstring, optopt )
    if ( i == 0 ) then
       ! unrecognized option
       process_short = '?'
       if ( opterr ) then
          print '(a,a,a)', "Error: unrecognized option '-", optopt, "'"
       endif
    endif
    if ( i > 0 .and. substr( optstring, i+1, i+1 ) == ':' ) then
       ! required argument
       optind = optind + 1
       if ( arglen > grpind ) then
          ! -xarg, return remainder of arg
          optarg = arg(grpind+1:arglen)
       elseif ( optind <= iargc()) then
          ! -x arg, return next arg
          call getarg( optind, optarg )
          optind = optind + 1
       elseif ( opterr ) then
          print '(a,a,a)', "Error: option '-", optopt, "' requires an argument"
       endif
       grpind = 2
    elseif ( arglen > grpind ) then
       ! no argument (or unrecognized), move to next option in argument (-xyz)
       grpind = grpind + 1
    else
       ! no argument (or unrecognized), move to to next argument
       grpind = 2
       optind = optind + 1
    endif
  end function process_short

end module getopt_m
