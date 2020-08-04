# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.1...3.15)
project( Optimal_bounding_box_Benchmark )

# CGAL and its components
find_package( CGAL QUIET )

if ( NOT CGAL_FOUND )
  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()
endif()

# include helper file
include( ${CGAL_USE_FILE} )

find_package(Eigen3 3.1.0 REQUIRED) #(3.1.0 or greater)
include(CGAL_Eigen_support)
if (NOT TARGET CGAL::Eigen_support)
  message(STATUS "This project requires the Eigen library, and will not be compiled.")
  return()
endif()

create_single_source_cgal_program("bench_obb.cpp")
target_link_libraries(bench_obb PUBLIC CGAL::Eigen_support)
