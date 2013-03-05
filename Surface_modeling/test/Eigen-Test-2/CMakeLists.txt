#cmake -D EIGEN3_DEVEL_INCLUDE=Path to Eigen dev level version
project(Eigen_SparseLU_Test)
cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
  cmake_policy(VERSION 2.8.4)
else()
  cmake_policy(VERSION 2.6)
endif()

# Find CGAL and CGAL Qt4
find_package(CGAL COMPONENTS Qt4)
include( ${CGAL_USE_FILE} )

# Find Eigen3 (requires 3.1.0 or greater)
find_package(Eigen3 3.1.0)


  if(EIGEN3_FOUND)
    include( ${EIGEN3_USE_FILE} )
  else()
    message(STATUS "NOTICE: Eigen library is not found.")
  endif()


#set(EIGEN3_DEVEL_INCLUDE "" CACHE FILEPATH "Path to Eigen dev level version")
#include_directories(${EIGEN3_DEVEL_INCLUDE})

#find_package(SuperLU)
#if (SUPERLU_FOUND)
#	include( ${SUPERLU_USE_FILE} )
#endif()

add_executable (Test main.cpp)
