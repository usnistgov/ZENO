# - Try to find sprng
# Once done this will define
#  SPRNG_FOUND - System has sprng
#  SPRNG_INCLUDE_DIRS - The sprng include directories
#  SPRNG_LIBRARIES - The libraries needed to use sprng
#  SPRNG_DEFINITIONS - Compiler switches required for using sprng

find_package(PkgConfig)
pkg_check_modules(PC_SPRNG QUIET sprng)
set(SPRNG_DEFINITIONS ${PC_SPRNG_CFLAGS_OTHER})

find_path(SPRNG_INCLUDE_DIR sprng.h
          HINTS ${PC_SPRNG_INCLUDEDIR} ${PC_SPRNG_INCLUDE_DIRS} )

find_library(SPRNG_LIBRARY NAMES sprng libsprng
             HINTS ${PC_SPRNG_LIBDIR} ${PC_SPRNG_LIBRARY_DIRS} )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set SPRNG_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(SPRNG  DEFAULT_MSG
                                  SPRNG_LIBRARY SPRNG_INCLUDE_DIR)

mark_as_advanced(SPRNG_INCLUDE_DIR SPRNG_LIBRARY )

set(SPRNG_LIBRARIES ${SPRNG_LIBRARY} )
set(SPRNG_INCLUDE_DIRS ${SPRNG_INCLUDE_DIR} )
