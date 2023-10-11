# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.1...3.23)
project(Optimal_bounding_box_Benchmark)

# CGAL and its components
find_package(CGAL REQUIRED)

find_package(Eigen3 3.1.0 QUIET) #(3.1.0 or greater)
include(CGAL_Eigen3_support)
if(NOT TARGET CGAL::Eigen3_support)
  message("NOTICE: This project requires the Eigen library, and will not be compiled.")
  return()
endif()

create_single_source_cgal_program("bench_obb.cpp")
target_link_libraries(bench_obb PUBLIC CGAL::Eigen3_support)
