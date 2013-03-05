#cmake -D EIGEN3_DEVEL_INCLUDE=Path to Eigen dev level version
project(Eigen_SparseLU_Test)
cmake_minimum_required(VERSION 2.6.2)

set(EIGEN3_DEVEL_INCLUDE "" CACHE FILEPATH "Path to Eigen dev level version")
include_directories(${EIGEN3_DEVEL_INCLUDE})

#find_package(SuperLU)
#if (SUPERLU_FOUND)
#	include( ${SUPERLU_USE_FILE} )
#endif()

add_executable (Test main.cpp)
