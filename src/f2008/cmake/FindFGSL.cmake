# ================================================================
# 
# Author: Walid Keyrouz <walid.keyrouz@nist.gov>
# Date:   Wed Jan 14 13:28:32 2015 EST
# 
# Time-stamp: <2015-01-15 15:06:16 wtk>
# 
# Notes:
# 
# - Try to find FGSL.  When done, this will define
# 
#   FGSL_FOUND        - system has FGSL (presumably also GSL)
#   FGSL_INCLUDE_DIRS - where to find fgsl.mod
#   FGSL_LIBRARIRES   - link these to use FGSL & GSL
# 
# - Based on http://www.cmake.org/Wiki/CMake:How_To_Find_Libraries
# 
# ================================================================

include (LibFindMacros)

# Dependencies
libfind_package (FGSL GSL REQUIRED)

# Use pkg-config to get hints about paths
libfind_pkg_check_modules (FGSL_PKGCONF fgsl)

# Include dir
find_path (FGSL_INCLUDE_DIR
  NAMES fgsl.mod
  PATHS ${FGSL_PKGCONF_INCLUDE_DIRS}
  )

# The library itself
find_library (FGSL_LIBRARY
  NAMES fgsl
  PATHS ${FGSL_PKGCONFIG_LIBRARY_DIRS}
  )

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
set (FGSL_PROCESS_INCLUDES FGSL_INCLUDE_DIR)
set (FGSL_PROCESS_LIBS FGSL_LIBRARY)

libfind_process (FGSL)

if (0)

  message ("FGSL_FOUND: " ${FGSL_FOUND})
  message ("FGSL_INCLUDE_DIRS: " "${FGSL_INCLUDE_DIRS}")
  message ("FGSL_LIBRARIES: " "${FGSL_LIBRARIES}")

endif (0)

# ================================================================

# Local Variables:
# mode: cmake-mode
# time-stamp-line-limit: 30
# End:
