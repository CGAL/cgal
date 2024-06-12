# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.12...3.29)
project(Surface_mesh_skeletonization_Tests)

find_package(CGAL REQUIRED)

find_package(Eigen3 3.2.0 QUIET) #(requires 3.2.0 or greater)
include(CGAL_Eigen3_support)
if(TARGET CGAL::Eigen3_support)
  create_single_source_cgal_program("MCF_Skeleton_test.cpp")
  target_link_libraries(MCF_Skeleton_test PUBLIC CGAL::Eigen3_support)
  create_single_source_cgal_program("skeleton_connectivity_test.cpp")
  target_link_libraries(skeleton_connectivity_test PUBLIC CGAL::Eigen3_support)
else()
  message("NOTICE: These tests require the Eigen library (3.2 or greater), and will not be compiled.")
endif()
