# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.1...3.23)
project(Optimal_bounding_box_Tests)

find_package(CGAL REQUIRED)

find_package(Eigen3 3.1.0 QUIET) #(3.1.0 or greater)
include(CGAL_Eigen3_support)
if(NOT TARGET CGAL::Eigen3_support)
  message("NOTICE: This project requires the Eigen library, and will not be compiled.")
  return()
endif()

create_single_source_cgal_program("test_OBB_traits.cpp")
create_single_source_cgal_program("test_nelder_mead.cpp")
create_single_source_cgal_program("test_optimization_algorithms.cpp")

foreach(target test_OBB_traits test_nelder_mead test_optimization_algorithms)
  target_link_libraries(${target} PUBLIC CGAL::Eigen3_support)
endforeach()
