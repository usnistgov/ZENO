# ================================================================
# 
# Author: Walid Keyrouz <walid.keyrouz@nist.gov>
# Date:   Wed Jan 14 13:28:32 2015 EST
# 
# Time-stamp: <2015-01-15 12:42:21 wtk>
# 
# Notes:
# 
# - Try to find GSL.  When done, this will define
# 
#   GSL_FOUND        - system has GSL
#   GSL_INCLUDE_DIRS - where to find gsl.mod
#   GSL_LIBRARIRES   - link these to use GSL
# 
# - Based on http://www.cmake.org/Wiki/CMake:How_To_Find_Libraries
# 
# ================================================================

include (LibFindMacros)

# Use pkg-config to get hints about paths
libfind_pkg_check_modules (GSL_PKGCONF gsl)

# Include dir
find_path (GSL_INCLUDE_DIR
  NAMES gsl/gsl_version.h
  PATHS ${GSL_PKGCONF_INCLUDE_DIRS}
  )

# The libraries itself
find_library (GSL_LIBRARY
  NAMES gsl gslcblas
  PATHS ${GSL_PKGCONF_LIBRARY_DIRS}
  )

find_library (GSL_CBLAS_LIBRARY
  NAMES gslcblas
  PATHS ${GSL_PKGCONF_LIBRARY_DIRS}
  )

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
set (GSL_PROCESS_INCLUDES GSL_INCLUDE_DIR)
set (GSL_PROCESS_LIBS GSL_LIBRARY GSL_CBLAS_LIBRARY)

libfind_process (GSL)

if (0)

  message ("GSL_FOUND: " ${GSL_FOUND})
  message ("GSL_INCLUDE_DIRS: " "${GSL_INCLUDE_DIRS}")
  message ("GSL_LIBRARIES: " "${GSL_LIBRARIES}")

endif (0)

# ================================================================

# Local Variables:
# mode: cmake-mode
# time-stamp-line-limit: 30
# End:
