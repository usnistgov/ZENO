cmake_minimum_required(VERSION 3.1)

set(CMAKE_BUILD_TYPE_INIT "Release")

project(ZENO VERSION 5.2.0)

include(GNUInstallDirs)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED on)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")

set(SRC_DIR ${CMAKE_SOURCE_DIR}/src)

set(CMDLINE_DIR ${SRC_DIR})
set(CMDLINE_C cmdline.c)
set(CMDLINE_H cmdline.h)

set(BOD_PARSER_DIR ${SRC_DIR}/BodParser)
set(BOD_PARSE_CC bod_parse.cc)
set(BOD_PARSERBASE_H BodParserbase.h)

set(BOD_SCANNER_DIR ${SRC_DIR}/BodParser)
set(BOD_SCAN_CC bod_lex.cc)
set(BOD_SCANNERBASE_H BodScannerbase.h)

set(SOURCE_FILES_LIB
  ${SRC_DIR}/Zeno.cc
  ${SRC_DIR}/ParametersWalkOnSpheres.cc
  ${SRC_DIR}/ParametersInteriorSampling.cc
  ${SRC_DIR}/ParametersResults.cc
  ${SRC_DIR}/Results.cc
  ${SRC_DIR}/ResultsCompiler.cc
  ${SRC_DIR}/ResultsInterior.cc
  ${SRC_DIR}/ResultsZeno.cc
  ${SRC_DIR}/Timer.cc
  ${SRC_DIR}/Units.cc)

set(HEADER_FILES_LIB
  ${SRC_DIR}/BoundingSphereGenerator.h
  ${SRC_DIR}/CovarianceMatrix.h
  ${SRC_DIR}/ParametersInteriorSampling.h
  ${SRC_DIR}/ParametersResults.h
  ${SRC_DIR}/ParametersWalkOnSpheres.h
  ${SRC_DIR}/Result.h
  ${SRC_DIR}/ResultsCompiler.h
  ${SRC_DIR}/Results.h
  ${SRC_DIR}/ResultsInterior.h
  ${SRC_DIR}/ResultsZeno.h
  ${SRC_DIR}/Timer.h
  ${SRC_DIR}/Uncertain.h
  ${SRC_DIR}/Units.h
  ${SRC_DIR}/Zeno.h
  ${SRC_DIR}/Geometry/AABBTree.h
  ${SRC_DIR}/Geometry/Cuboid.h
  ${SRC_DIR}/Geometry/Matrix3x3.h
  ${SRC_DIR}/Geometry/MixedModel.h
  ${SRC_DIR}/Geometry/MixedModelProcessed.h
  ${SRC_DIR}/Geometry/Sphere.h
  ${SRC_DIR}/Geometry/Triangle.h
  ${SRC_DIR}/Geometry/Vector3.h
  ${SRC_DIR}/Geometry/Voxels.h
  ${SRC_DIR}/RandomNumber/LeapFrog.h
  ${SRC_DIR}/RandomNumber/Rand.h
  ${SRC_DIR}/RandomNumber/SPRNG.h
  ${SRC_DIR}/SpherePoint/BiasedSpherePointDirect.h
  ${SRC_DIR}/SpherePoint/BiasedSpherePointRejection.h
  ${SRC_DIR}/SpherePoint/QuasiRandomSpherePointsRegularCube.h
  ${SRC_DIR}/SpherePoint/RandomBallPointRejection.h
  ${SRC_DIR}/SpherePoint/RandomSpherePointMarsaglia.h
  ${SRC_DIR}/SpherePoint/RandomSpherePointPolar.h
  ${SRC_DIR}/Walker/SamplerInterior.h
  ${SRC_DIR}/Walker/WalkerExterior.h)

set(SOURCE_FILES_CLI
  ${CMDLINE_DIR}/${CMDLINE_C}
  ${BOD_PARSER_DIR}/${BOD_PARSE_CC}
  ${BOD_PARSER_DIR}/BodParser.cc
  ${BOD_SCANNER_DIR}/${BOD_SCAN_CC}
  ${SRC_DIR}/Main.cc
  ${SRC_DIR}/ParametersLocal.cc)

set(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})

add_library(libzeno SHARED ${SOURCE_FILES_LIB})

set_target_properties(libzeno PROPERTIES
  VERSION ${PROJECT_VERSION}
  PUBLIC_HEADER "${HEADER_FILES_LIB}"
  OUTPUT_NAME zeno)

install(TARGETS libzeno
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

add_executable(zeno-cli ${SOURCE_FILES_CLI})

set_target_properties(zeno-cli PROPERTIES
  OUTPUT_NAME zeno)

target_link_libraries(zeno-cli LINK_PUBLIC libzeno)

install(TARGETS zeno-cli
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

find_package(ZLIB REQUIRED)
target_link_libraries(libzeno ${ZLIB_LIBRARIES})

find_package(Threads REQUIRED)
target_link_libraries(libzeno ${CMAKE_THREAD_LIBS_INIT})

#setup rng

option(USE_SPRNG "Use SPRNG random number generator" ON)
mark_as_advanced(USE_SPRNG)

if (USE_SPRNG)
  find_package(SPRNG)

  if (NOT SPRNG_FOUND)
    mark_as_advanced(CLEAR SPRNG_INCLUDE_DIR)
    mark_as_advanced(CLEAR SPRNG_LIBRARY)

    message(SEND_ERROR "SPRNG not found!")
  else()
    include_directories(${SPRNG_INCLUDE_DIRS})
    target_link_libraries(libzeno ${SPRNG_LIBRARIES})

    add_definitions(-DUSE_SPRNG_RNG)

    target_sources(libzeno PUBLIC ${SRC_DIR}/RandomNumber/SPRNG.cc)
  endif()
else()
  message(WARNING "USE_SPRNG is OFF.  Results have not been validated without the use of SPRNG!")
  
  add_definitions(-DUSE_LEAP_FROG_RNG)

  target_sources(libzeno PUBLIC ${SRC_DIR}/RandomNumber/LeapFrog.cc)
endif()

#setup mpi

option(USE_MPI "Use MPI" OFF)

if (USE_MPI)
  find_package(MPI REQUIRED)
  
  include_directories(${MPI_INCLUDE_PATH})
  target_link_libraries(libzeno ${MPI_LIBRARIES})

  add_definitions(-DUSE_MPI)
  
  if(MPI_COMPILE_FLAGS)
    set_target_properties(libzeno PROPERTIES
      COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
  endif()

  if(MPI_LINK_FLAGS)
    set_target_properties(libzeno PROPERTIES
      LINK_FLAGS "${MPI_LINK_FLAGS}")
  endif()
endif()

#setup gengetopt

option (GEN_CMD_PARSER "Generate command line parser" OFF)
mark_as_advanced(GEN_CMD_PARSER)

if (GEN_CMD_PARSER)
  find_program(GENGETOPT_BIN gengetopt)

  mark_as_advanced(GENGETOPT_BIN)

  if (GENGETOPT_BIN)
    set(OPTIONS_GGO ${CMDLINE_DIR}/options.ggo)
    
    add_custom_command(
      OUTPUT
      ${CMDLINE_DIR}/${CMDLINE_C}
      ${CMDLINE_DIR}/${CMDLINE_H}
      COMMAND
      ${GENGETOPT_BIN}
      -i ${OPTIONS_GGO}
      --conf-parser
      --output-dir=${CMDLINE_DIR}
      --set-version=${PROJECT_VERSION}
      DEPENDS ${OPTIONS_GGO}
      VERBATIM)
  else()
    message(SEND_ERROR "gengetopt not found!")
  endif()
endif()

#setup file parsers

macro(setup_file_parser
    OPTION_NAME
    OPTION_DESCRIPTION
    NAMESPACE
    PARSER_DIR
    PARSE_CC
    PARSERBASE_H
    PARSER_CLASS
    GRAMMAR
    SCANNER_DIR
    SCAN_CC
    SCANNERBASE_H
    SCANNER_CLASS
    LEXER)
  
  option(${OPTION_NAME} ${OPTION_DESCRIPTION} OFF)
  mark_as_advanced(${OPTION_NAME})

  if (${OPTION_NAME})
    #setup bisonc++
    
    find_program(BISONCXX_BIN bisonc++)

    mark_as_advanced(BISONCXX_BIN)

    if (BISONCXX_BIN)
      add_custom_command(
	OUTPUT ${PARSER_DIR}/${PARSE_CC} ${PARSER_DIR}/${PARSERBASE_H}
	COMMAND ${BISONCXX_BIN}
	--target-directory ${PARSER_DIR}
	--class-name ${PARSER_CLASS}
	--scanner-class-name ${SCANNER_CLASS}
	--parsefun-source ${PARSE_CC}
	--namespace ${NAMESPACE}
	${PARSER_DIR}/${GRAMMAR}
	DEPENDS ${PARSER_DIR}/${GRAMMAR}
	VERBATIM)
    else()
      message(SEND_ERROR "bisonc++ not found!")
    endif()

    #setup flexc++
    
    find_program(FLEXCXX_BIN flexc++)

    mark_as_advanced(FLEXCXX_BIN)

    if (FLEXCXX_BIN)
      add_custom_command(
	OUTPUT ${SCANNER_DIR}/${SCAN_CC} ${SCANNER_DIR}/${SCANNERBASE_H}
	COMMAND ${FLEXCXX_BIN}
	--target-directory ${SCANNER_DIR}
	--class-name ${SCANNER_CLASS}
	--lex-source ${SCAN_CC}
	--namespace ${NAMESPACE}
	${SCANNER_DIR}/${LEXER}
	DEPENDS ${SCANNER_DIR}/${LEXER}
	VERBATIM)
    else()
      message(SEND_ERROR "flexc++ not found!")
    endif()
  endif()
endmacro(setup_file_parser)

setup_file_parser(
  GEN_BOD_PARSER
  "Generate bod file parser"
  bod_parser
  ${BOD_PARSER_DIR}
  ${BOD_PARSE_CC}
  ${BOD_PARSERBASE_H}
  BodParser
  bod_grammar
  ${BOD_SCANNER_DIR}
  ${BOD_SCAN_CC}
  ${BOD_SCANNERBASE_H}
  BodScanner
  bod_lexer)

#setup testing

set(TEST_SCRIPT ${SRC_DIR}/SelfTests/run_tests.py)
set(TEST_EXECUTABLE ${CMAKE_BINARY_DIR}/zeno)
set(TEST_INPUT_DIR ${SRC_DIR}/SelfTests/Input)
set(TEST_GROUND_TRUTH_DIR ${SRC_DIR}/SelfTests/GroundTruth)
set(TEST_OUTPUT_DIR ${CMAKE_BINARY_DIR}/SelfTestOutput)

if (USE_MPI)
  set(TEST_FLAGS --run-mpi-tests)
else()
  set(TEST_FLAGS)
endif()

file(MAKE_DIRECTORY ${TEST_OUTPUT_DIR})

add_custom_target(check 
                  COMMAND ${TEST_SCRIPT} ${TEST_FLAGS} ${TEST_EXECUTABLE} ${TEST_INPUT_DIR} ${TEST_GROUND_TRUTH_DIR} ${TEST_OUTPUT_DIR}
                  VERBATIM)

