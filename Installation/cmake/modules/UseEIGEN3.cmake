# This module setups the compiler for using EIGEN3 library.
# It assumes that find_package(EIGEN3) was already called.


include_directories( ${EIGEN3_INCLUDE_DIR} )
include_directories( ${EIGEN3_INCLUDE_DIR}/unsupported )

add_definitions(-DCGAL_EIGEN3_ENABLED)
