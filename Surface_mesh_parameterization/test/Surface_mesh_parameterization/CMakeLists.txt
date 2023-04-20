# This is the CMake script for compiling this folder.

cmake_minimum_required(VERSION 3.1...3.23)
project(Surface_mesh_parameterization_Tests)

# Find CGAL
find_package(CGAL REQUIRED)

find_package(Eigen3 3.1.0 QUIET) #(requires 3.1.0 or greater)
include(CGAL_Eigen3_support)
if(TARGET CGAL::Eigen3_support)
  create_single_source_cgal_program("extensive_parameterization_test.cpp")
  target_link_libraries(extensive_parameterization_test PUBLIC CGAL::Eigen3_support)
else()
  message("NOTICE: The tests require Eigen 3.1 (or greater), and will not be compiled.")
endif()
