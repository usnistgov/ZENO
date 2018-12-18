cmake_minimum_required(VERSION 3.1)

set(CMAKE_BUILD_TYPE_INIT "Release")

project(Zeno)

#set(CMAKE_VERBOSE_MAKEFILE on)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED on)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")

set(SRC_DIR ${CMAKE_SOURCE_DIR}/src)

set(CMDLINE_DIR ${SRC_DIR})
set(CMDLINE_C ${CMDLINE_DIR}/cmdline.c)
set(CMDLINE_H ${CMDLINE_DIR}/cmdline.h)

set(PARSER_DIR ${SRC_DIR})
set(PARSE_CC ${PARSER_DIR}/parse.cc)
set(PARSERBASE_H ${PARSER_DIR}/Parserbase.h)

set(SCANNER_DIR ${SRC_DIR})
set(SCAN_CC ${SCANNER_DIR}/lex.cc)
set(SCANNERBASE_H ${SCANNER_DIR}/Scannerbase.h)

set(SOURCE_FILES
  ${CMDLINE_C}
  ${PARSE_CC}
  ${SCAN_CC}
  ${SRC_DIR}/Main.cc
  ${SRC_DIR}/Parameters.cc
  ${SRC_DIR}/Parser.cc
  ${SRC_DIR}/ResultsCompiler.cc
  ${SRC_DIR}/ResultsInterior.cc
  ${SRC_DIR}/ResultsZeno.cc
  ${SRC_DIR}/Timer.cc
  ${SRC_DIR}/Units.cc)

add_executable(zeno ${SOURCE_FILES})

install(TARGETS zeno
        RUNTIME DESTINATION bin)

find_package(ZLIB REQUIRED)
target_link_libraries(zeno ${ZLIB_LIBRARIES})

find_package(Threads REQUIRED)
target_link_libraries(zeno ${CMAKE_THREAD_LIBS_INIT})

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
    target_link_libraries(zeno ${SPRNG_LIBRARIES})

    add_definitions(-DUSE_SPRNG_RNG)

    target_sources(zeno PUBLIC ${SRC_DIR}/RandomNumber/SPRNG.cc)
  endif()
else()
  add_definitions(-DUSE_LEAP_FROG_RNG)

  target_sources(zeno PUBLIC ${SRC_DIR}/RandomNumber/LeapFrog.cc)
endif()

#setup mpi

option(USE_MPI "Use MPI" OFF)

if (USE_MPI)
  find_package(MPI REQUIRED)
  
  include_directories(${MPI_INCLUDE_PATH})
  target_link_libraries(zeno ${MPI_LIBRARIES})

  add_definitions(-DUSE_MPI)
  
  if(MPI_COMPILE_FLAGS)
    set_target_properties(zeno PROPERTIES
      COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
  endif()

  if(MPI_LINK_FLAGS)
    set_target_properties(zeno PROPERTIES
      LINK_FLAGS "${MPI_LINK_FLAGS}")
  endif()
endif()

#setup gengetopt

find_program(GENGETOPT_BIN gengetopt)

mark_as_advanced(GENGETOPT_BIN)

if (GENGETOPT_BIN)
  set(OPTIONS_GGO ${CMDLINE_DIR}/options.ggo)
  
  add_custom_command(OUTPUT  ${CMDLINE_C} ${CMDLINE_H}
                     COMMAND ${GENGETOPT_BIN} -i ${OPTIONS_GGO} --conf-parser --output-dir=${CMDLINE_DIR}
                     DEPENDS ${OPTIONS_GGO}
                     VERBATIM)
endif()

#setup bisonc++

find_program(BISONCXX_BIN bisonc++)

mark_as_advanced(BISONCXX_BIN)

if (BISONCXX_BIN)
  set(GRAMMAR ${PARSER_DIR}/grammar)
  
  add_custom_command(OUTPUT  ${PARSE_CC} ${PARSERBASE_H}
                     COMMAND ${BISONCXX_BIN} --target-directory ${PARSER_DIR} ${GRAMMAR}
                     DEPENDS ${GRAMMAR}
                     VERBATIM)
endif()

#setup flexc++

find_program(FLEXCXX_BIN flexc++)

mark_as_advanced(FLEXCXX_BIN)

if (FLEXCXX_BIN)
  set(LEXER ${SCANNER_DIR}/lexer)
  
  add_custom_command(OUTPUT  ${SCAN_CC} ${SCANNERBASE_H}
                     COMMAND ${FLEXCXX_BIN} --target-directory ${SCANNER_DIR} ${LEXER}
                     DEPENDS ${LEXER}
                     VERBATIM)
endif()

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

