# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( optimal_rotation_benchmark)

cmake_minimum_required(VERSION 2.6.2)

find_package(CGAL QUIET COMPONENTS Core )


if ( CGAL_FOUND )

  add_definitions(-DDEBUG_TRACE)

  include( ${CGAL_USE_FILE} )

  add_executable  ( optimal_rotation_svd_eigen
	optimal_rotation_svd_eigen.cpp)

  add_executable  ( optimal_rotation_svd_nr
	optimal_rotation_svd_nr.cpp
	nrutil.cpp)

  add_executable  ( optimal_rotation_polar
	optimal_rotation_polar.cpp)

  add_executable  ( optimal_rotation_polar_eigen
	optimal_rotation_polar_eigen.cpp)

else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()


find_package(Eigen3 3.0.0)

if(EIGEN3_FOUND)
    include( ${EIGEN3_USE_FILE} )
  endif(EIGEN3_FOUND)



