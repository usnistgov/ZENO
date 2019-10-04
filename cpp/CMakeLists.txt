cmake_minimum_required(VERSION 3.1)

set(CMAKE_BUILD_TYPE_INIT "Release")

project(ZENO VERSION 5.2.1)

include(GNUInstallDirs)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED on)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")

set(SRC_DIR_LIB ${CMAKE_SOURCE_DIR}/libzeno/src)
set(INCLUDE_DIR_LIB ${CMAKE_SOURCE_DIR}/libzeno/include)

set(SRC_DIR_CLI ${CMAKE_SOURCE_DIR}/zeno-cli/src)

set(CMDLINE_DIR ${SRC_DIR_CLI})
set(CMDLINE_C cmdline.c)
set(CMDLINE_H cmdline.h)

set(BOD_PARSER_DIR ${SRC_DIR_CLI}/BodParser)
set(BOD_PARSE_CC bod_parse.cc)
set(BOD_PARSERBASE_H BodParserbase.h)

set(BOD_SCANNER_DIR ${SRC_DIR_CLI}/BodParser)
set(BOD_SCAN_CC bod_lex.cc)
set(BOD_SCANNERBASE_H BodScannerbase.h)

set(SOURCE_FILES_LIB
  ${SRC_DIR_LIB}/Zeno.cc
  ${SRC_DIR_LIB}/ParametersWalkOnSpheres.cc
  ${SRC_DIR_LIB}/ParametersInteriorSampling.cc
  ${SRC_DIR_LIB}/ParametersResults.cc
  ${SRC_DIR_LIB}/Results.cc
  ${SRC_DIR_LIB}/ResultsCompiler.cc
  ${SRC_DIR_LIB}/ResultsInterior.cc
  ${SRC_DIR_LIB}/ResultsZeno.cc
  ${SRC_DIR_LIB}/Timer.cc
  ${SRC_DIR_LIB}/Units.cc)

set(HEADER_FILES_LIB
  ${INCLUDE_DIR_LIB}/BoundingSphereGenerator.h
  ${INCLUDE_DIR_LIB}/CovarianceMatrix.h
  ${INCLUDE_DIR_LIB}/ParametersInteriorSampling.h
  ${INCLUDE_DIR_LIB}/ParametersResults.h
  ${INCLUDE_DIR_LIB}/ParametersWalkOnSpheres.h
  ${INCLUDE_DIR_LIB}/Result.h
  ${INCLUDE_DIR_LIB}/ResultsCompiler.h
  ${INCLUDE_DIR_LIB}/Results.h
  ${INCLUDE_DIR_LIB}/ResultsInterior.h
  ${INCLUDE_DIR_LIB}/ResultsZeno.h
  ${INCLUDE_DIR_LIB}/Timer.h
  ${INCLUDE_DIR_LIB}/Uncertain.h
  ${INCLUDE_DIR_LIB}/Units.h
  ${INCLUDE_DIR_LIB}/Zeno.h
  ${INCLUDE_DIR_LIB}/Geometry/AABBTree.h
  ${INCLUDE_DIR_LIB}/Geometry/Cuboid.h
  ${INCLUDE_DIR_LIB}/Geometry/Matrix3x3.h
  ${INCLUDE_DIR_LIB}/Geometry/MixedModel.h
  ${INCLUDE_DIR_LIB}/Geometry/MixedModelProcessed.h
  ${INCLUDE_DIR_LIB}/Geometry/Sphere.h
  ${INCLUDE_DIR_LIB}/Geometry/Triangle.h
  ${INCLUDE_DIR_LIB}/Geometry/Vector3.h
  ${INCLUDE_DIR_LIB}/Geometry/Voxels.h
  ${INCLUDE_DIR_LIB}/RandomNumber/LeapFrog.h
  ${INCLUDE_DIR_LIB}/RandomNumber/Rand.h
  ${INCLUDE_DIR_LIB}/RandomNumber/SPRNG.h
  ${INCLUDE_DIR_LIB}/SpherePoint/BiasedSpherePointDirect.h
  ${INCLUDE_DIR_LIB}/SpherePoint/BiasedSpherePointRejection.h
  ${INCLUDE_DIR_LIB}/SpherePoint/QuasiRandomSpherePointsRegularCube.h
  ${INCLUDE_DIR_LIB}/SpherePoint/RandomBallPointRejection.h
  ${INCLUDE_DIR_LIB}/SpherePoint/RandomSpherePointMarsaglia.h
  ${INCLUDE_DIR_LIB}/SpherePoint/RandomSpherePointPolar.h
  ${INCLUDE_DIR_LIB}/Walker/SamplerInterior.h
  ${INCLUDE_DIR_LIB}/Walker/WalkerExterior.h)

set(SOURCE_FILES_CLI
  ${CMDLINE_DIR}/${CMDLINE_C}
  ${BOD_PARSER_DIR}/${BOD_PARSE_CC}
  ${BOD_PARSER_DIR}/BodParser.cc
  ${BOD_SCANNER_DIR}/${BOD_SCAN_CC}
  ${SRC_DIR_CLI}/Main.cc
  ${SRC_DIR_CLI}/ParametersLocal.cc)

set(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})

add_library(libzeno SHARED ${SOURCE_FILES_LIB})

target_include_directories(libzeno PUBLIC ${INCLUDE_DIR_LIB})

set_target_properties(libzeno PROPERTIES
  VERSION ${PROJECT_VERSION}
  OUTPUT_NAME zeno)

install(TARGETS libzeno
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

install(DIRECTORY ${INCLUDE_DIR_LIB}/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

add_executable(zeno-cli ${SOURCE_FILES_CLI})

target_include_directories(zeno-cli PRIVATE ${INCLUDE_DIR_LIB})

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

    target_sources(libzeno PUBLIC ${SRC_DIR_LIB}/RandomNumber/SPRNG.cc)
  endif()
else()
  message(WARNING "USE_SPRNG is OFF.  Results have not been validated without the use of SPRNG!")
  
  add_definitions(-DUSE_LEAP_FROG_RNG)

  target_sources(libzeno PUBLIC ${SRC_DIR_LIB}/RandomNumber/LeapFrog.cc)
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

set(TEST_SCRIPT ${SRC_DIR_CLI}/SelfTests/run_tests.py)
set(TEST_EXECUTABLE ${CMAKE_BINARY_DIR}/zeno)
set(TEST_INPUT_DIR ${SRC_DIR_CLI}/SelfTests/Input)
set(TEST_GROUND_TRUTH_DIR ${SRC_DIR_CLI}/SelfTests/GroundTruth)
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

