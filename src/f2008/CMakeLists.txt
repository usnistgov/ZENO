# ================================================================
# 
# CMake project for zeno-gsl
# 
# Author: Walid Keyrouz <walid.keyrouz@nist.gov>
# Date: Fri Jan  9 15:40:18 2015 EST
#
# Time-stamp: <2015-08-31 15:42:35 walid>
# 
# ================================================================

cmake_minimum_required (VERSION 2.8)

# ----------------------------------------------------------------

# These options should appear before the PROJECT keyword
set (CMAKE_DISABLE_SOURCE_CHANGES ON)
set (CMAKE_DISABLE_IN_SOURCE_BUILD ON)

# ---------------------------------------------------------------- 

project (zeno-v4)
enable_language (Fortran)

# Enable private CMake modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

# ----------------------------------------------------------------

if (ZENO_CMAKE_DBG GREATER 3)
  message ("")
  message ("* Found ZENO_CMAKE_DBG")
endif (ZENO_CMAKE_DBG GREATER 3)

if (ZENO_CMAKE_DBG GREATER 2)
  message ("")
  message ("* Default CMAKE_Fortran_FLAGS_<build_type>")
  message ("* CMAKE_Fortran_FLAGS_DEBUG:   ${CMAKE_Fortran_FLAGS_DEBUG}")
  message ("* CMAKE_Fortran_FLAGS_RELEASE: ${CMAKE_Fortran_FLAGS_RELEASE}")
  message ("* CMAKE_Fortran_FLAGS_PROFILE: ${CMAKE_Fortran_FLAGS_PROFILE}")
  message ("* CMAKE_Fortran_FLAGS_RELWITHDEBINFO: "
    "${CMAKE_Fortran_FLAGS_RELWITHDEBINFO}")
endif (ZENO_CMAKE_DBG GREATER 2)

# ----------------------------------------------------------------

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
    "Choose the build type, one of: "
    "None Debug Release RelWithDebInfo Profile."
    FORCE)
endif (NOT CMAKE_BUILD_TYPE)

if (ZENO_CMAKE_DBG GREATER 1)
  message ("")
  message ("* Current build type: ${CMAKE_BUILD_TYPE}")
endif (ZENO_CMAKE_DBG GREATER 1)

# ----------------------------------------------------------------

# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

set (GFortran_Common_Flags
  "-cpp"
  "-Wno-compare-reals"
  "-fno-f2c"
  )

set (GFortran_Common_Flags_Dbg
  ""
  )

set (GFortran_Common_Flags_Rel
  "-DNDEBUG"
  "-O3"
  "-funroll-all-loops"
  "-ffast-math"
  )

set (GFortran_Xtras_DEBUG
  "-O0"
  "-DZENO_ENABLE_LIB_OPTS"
  "-DZENO_ENABLE_RNG_LOG"
  "-Wall" "-Wextra" "-Wsurprising" "-Wimplicit-interface"
  "-fmax-errors=1" "-fcheck=all" "-fbacktrace"
  "${GFortran_Common_Flags}"
  "${GFortran_Common_Flags_Dbg}"
  )

set (GFortran_Xtras_RELEASE
  "${GFortran_Common_Flags}"
  )

set (GFortran_Xtras_RELWITHDEBINFO
  "${GFortran_Common_Flags}"
  "${GFortran_Common_Flags_Dbg}"
  "-O3"
  )

set (GFortran_Xtras_PROFILE
  "${GFortran_Common_Flags}"
  "${GFortran_Common_Flags_Dbg}"
  "-g"
  "-pg"
  )

# Convert lists to strings
string(REPLACE ";" " " GFortran_Xtras_DEBUG   "${GFortran_Xtras_DEBUG}")
string(REPLACE ";" " " GFortran_Xtras_RELEASE "${GFortran_Xtras_RELEASE}")
string(REPLACE ";" " " GFortran_Xtras_PROFILE "${GFortran_Xtras_PROFILE}")
string(REPLACE ";" " " GFortran_Xtras_RELWITHDEBINFO
  "${GFortran_Xtras_RELWITHDEBINFO}")

if (ZENO_CMAKE_DBG GREATER 3)
  message ("")
  message ("* GFortran_Xtras_DEBUG:   ${GFortran_Xtras_DEBUG}")
  message ("* GFortran_Xtras_RELEASE: ${GFortran_Xtras_RELEASE}")
  message ("* GFortran_Xtras_PROFILE: ${GFortran_Xtras_PROFILE}")
  message ("* GFortran_Xtras_RELWITHDEBINFO: "
    "${GFortran_Xtras_RELWITHDEBINFO}")

endif (ZENO_CMAKE_DBG GREATER 3)

# ----------------------------------------------------------------

if (Fortran_COMPILER_NAME MATCHES "gfortran.*")

  # gfortran
  set (CMAKE_Fortran_FLAGS_DEBUG
    "${CMAKE_Fortran_FLAGS_DEBUG} ${GFortran_Xtras_DEBUG}")
  set (CMAKE_Fortran_FLAGS_RELEASE
    "${CMAKE_Fortran_FLAGS_RELEASE} ${GFortran_Xtras_RELEASE}")
  set (CMAKE_Fortran_FLAGS_RELWITHDEBINFO
    "${CMAKE_Fortran_FLAGS_RELWITHDEBINFO} ${GFortran_Xtras_RELWITHDEBINFO}")
  set (CMAKE_Fortran_FLAGS_PROFILE
    "${CMAKE_Fortran_FLAGS_PROFILE} ${GFortran_Xtras_PROFILE}")

elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")

  # ifort (untested)
  set (CMAKE_Fortran_FLAGS_DEBUG   "-f77rtl -O0 -g")
  set (CMAKE_Fortran_FLAGS_RELEASE "-f77rtl -O3")

elseif (Fortran_COMPILER_NAME MATCHES "g77")

  # g77
  set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g -m32")
  set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3 -m32")

else ()

  message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
  message ("No optimized Fortran compiler flags are known, we just try -O2...")

  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set (CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-O2 -g")

endif (Fortran_COMPILER_NAME MATCHES "gfortran.*")

set (CMAKE_EXE_LINKER_FLAGS_PROFILE "-pg"
  CACHE STRING "Flags used for linking libraries during profile builds."
  FORCE
  )

set (CMAKE_SHARED_LINKER_FLAGS_PROFILE "-pg"
  CACHE STRING "Flags used by the shared libraries linker during profile builds."
  FORCE
  )

MARK_AS_ADVANCED (
  CMAKE_Fortran_FLAGS_PROFILE
  CMAKE_EXE_LINKER_FLAGS_PROFILE
  CMAKE_SHARED_LINKER_FLAGS_PROFILE
)

if (ZENO_CMAKE_DBG GREATER 1)
  message ("")
  message ("* Fortran_COMPILER_NAME:       ${Fortran_COMPILER_NAME}")
  message ("* CMAKE_Fortran_FLAGS_DEBUG:   ${CMAKE_Fortran_FLAGS_DEBUG}")
  message ("* CMAKE_Fortran_FLAGS_RELEASE: ${CMAKE_Fortran_FLAGS_RELEASE}")
  message ("* CMAKE_Fortran_FLAGS_PROFILE: ${CMAKE_Fortran_FLAGS_PROFILE}")
  message ("* CMAKE_Fortran_FLAGS_RELWITHDEBINFO: "
    "${CMAKE_Fortran_FLAGS_RELWITHDEBINFO}")
  message ("")
  message ("* CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
  message ("")
endif (ZENO_CMAKE_DBG GREATER 1)

# ================================================================

set (zeno_SRCS 
  getopt.f90 
  math-utils.f90 
  numeric-constants.f90 
  numeric-kinds.f90 
  sort-utils.f90 
  stat-utils.f90 
  zeno-codes-data.f90 
  zeno-debug.f90 
  zeno-eigen.f90 
  zeno-ellipsoid-utils.f90 
  zeno-enums.f90 
  zeno-files-data.f90 
  zeno-integrations.f90 
  zeno-kernels.f90 
  zeno-main.f90 
  zeno-options.f90 
  zeno-output.f90 
  zeno-parser.f90 
  zeno-prob-geom.f90 
  zeno-rng.f90 
  zeno-setup.f90 
  zeno-shape-dist.f90 
  zeno-shape-interior-test.f90 
  zeno-shape-surf-area.f90 
  zeno-shape-surf-points.f90 
  zeno-shape-volume.f90 
  zeno-sphere.f90 
  zeno-tensors.f90 
  zeno-toss.f90 
  zeno-units.f90 
  zeno-vectors.f90 
  zeno-version.f90
  zeno-warnings.f90 
  zeno-zeerot.f90 )

# ----------------------------------------------------------------

add_executable (zeno-v4 ${zeno_SRCS})

# ----------------------------------------------------------------

find_package (FGSL REQUIRED)
include_directories (${FGSL_INCLUDE_DIRS})
target_link_libraries (zeno-v4 ${FGSL_LIBRARIES})

# ================================================================

# Local Vars:
# mode: cmake-mode
# time-stamp-line-limit: 30
# End:
