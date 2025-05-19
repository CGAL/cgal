# This file sets up OSQP for CMake. Once done this will define
#  OSQP_FOUND       - system has OSQP lib
#  OSQP_INCLUDE_DIR - the OSQP include directory
#  OSQP_LIBRARIES   - link these to use OSQP

if(NOT OSQP_FOUND)
  find_path(OSQP_INCLUDE_DIR
    NAMES osqp.h
    PATHS /usr/local/include/osqp/
    ENV OSQP_INC_DIR)

  find_library(OSQP_LIBRARIES
    NAMES libosqp osqp
    PATHS ENV LD_LIBRARY_PATH
    ENV LIBRARY_PATH
    /usr/local/lib
    ${OSQP_INCLUDE_DIR}/../lib
    ENV OSQP_LIB_DIR)

  if(OSQP_LIBRARIES AND OSQP_INCLUDE_DIR)
    set(OSQP_FOUND TRUE)
  endif()
endif()
