# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.1...3.15)
project( Optimal_bounding_box_Tests )

find_package(CGAL QUIET)

if (NOT CGAL_FOUND)
  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()
endif()

include( ${CGAL_USE_FILE} )

find_package(Eigen3 3.1.0 REQUIRED) #(3.1.0 or greater)
if (NOT EIGEN3_FOUND)
  message(STATUS "This project requires the Eigen library, and will not be compiled.")
  return()
else()
  include(${EIGEN3_USE_FILE})
endif()

create_single_source_cgal_program("test_linear_algebra_functions.cpp")
create_single_source_cgal_program("test_optimization_algorithms.cpp")

