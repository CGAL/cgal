# - Check if the source file provided in the FILE argument compiles and runs.
# CHECK_CXX_FILE_RUNS(FILE VAR TEST)
#  FILE     - source file to try to compile
#  VAR      - variable set to 0 on success
#  TEST     - name of the test for a human being
#
# The following variables may be set before calling this macro to
# modify the way the check is run:
#
#  CMAKE_REQUIRED_FLAGS = string of compile command line flags
#  CMAKE_REQUIRED_DEFINITIONS = list of macros to define (-DFOO=bar)
#  CMAKE_REQUIRED_INCLUDES = list of include directories
#  CMAKE_REQUIRED_LIBRARIES = list of libraries to link
#
# Adapted to CGAL from IST's CheckCXXSourceRuns.cmake and
# KDE4's CheckCSourceRuns.cmake

MACRO(CHECK_CXX_FILE_RUNS FILE VAR TEST)

  # Set compiler settings
  SET(MACRO_CHECK_FUNCTION_DEFINITIONS "${CMAKE_REQUIRED_FLAGS}")
  if(CMAKE_REQUIRED_LIBRARIES)
    SET(CHECK_CXX_SOURCE_COMPILES_ADD_LIBRARIES "-DLINK_LIBRARIES:STRING=${CMAKE_REQUIRED_LIBRARIES}")
  else()
    SET(CHECK_CXX_SOURCE_COMPILES_ADD_LIBRARIES)
  endif()
  
  if(CMAKE_REQUIRED_INCLUDES)
    SET(CHECK_CXX_SOURCE_COMPILES_ADD_INCLUDES "-DINCLUDE_DIRECTORIES:STRING=${CMAKE_REQUIRED_INCLUDES}")
  else()
    SET(CHECK_CXX_SOURCE_COMPILES_ADD_INCLUDES)
  endif()

  # Try to compile and run the test
  #MESSAGE(STATUS "Performing Test ${TEST}")
  TRY_RUN(${VAR} ${VAR}_COMPILED
          "${CMAKE_BINARY_DIR}"
          "${FILE}"
          COMPILE_DEFINITIONS ${CMAKE_REQUIRED_DEFINITIONS}
          CMAKE_FLAGS -DCOMPILE_DEFINITIONS:STRING=${MACRO_CHECK_FUNCTION_DEFINITIONS}
          "${CHECK_CXX_SOURCE_COMPILES_ADD_LIBRARIES}"
          "${CHECK_CXX_SOURCE_COMPILES_ADD_INCLUDES}"
          OUTPUT_VARIABLE OUTPUT)
 
  # if it did not compile make the return value fail code of 1
  if(NOT ${VAR}_COMPILED)
    SET(${VAR} 1 CACHE INTERNAL "Test ${TEST}" FORCE )
  endif()
  # if the return value was 0 then it worked
  SET(result_var ${${VAR}})
  if("${result_var}" EQUAL 0)
    SET(${VAR} 1 CACHE INTERNAL "Test ${TEST}" FORCE )
    MESSAGE(STATUS "Performing Test ${TEST} - Success")
    FILE(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeOutput.log
      "Performing C++ SOURCE FILE Test ${TEST} succeeded with the following output:\n"
      "${OUTPUT}\n"
      "Source file was:\n${SOURCE}\n")
  else()
    MESSAGE(STATUS "Performing Test ${TEST} - Failed")
    SET(${VAR} "" CACHE INTERNAL "Test ${TEST}" FORCE )
    FILE(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log
      "Performing C++ SOURCE FILE Test ${TEST} failed with the following output:\n"
      "${OUTPUT}\n"
      "Return value: ${result_var}\n"
      "Source file was:\n${SOURCE}\n")
  endif()
ENDMACRO(CHECK_CXX_FILE_RUNS)

