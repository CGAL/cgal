# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.12...3.29)
project(Periodic_3_mesh_3_Tests)

find_package(CGAL REQUIRED COMPONENTS ImageIO)

# Use Eigen
find_package(Eigen3 3.1.0 QUIET) #(3.1.0 or greater)
include(CGAL_Eigen3_support)
if(NOT TARGET CGAL::Eigen3_support)
  message("NOTICE: This project requires the Eigen library, and will not be compiled.")
  return()
endif()

create_single_source_cgal_program("test_implicit_shapes_bunch.cpp")
target_link_libraries(test_implicit_shapes_bunch PUBLIC CGAL::Eigen3_support)
create_single_source_cgal_program("test_implicit_shapes_with_features.cpp")
target_link_libraries(test_implicit_shapes_with_features PUBLIC CGAL::Eigen3_support)
create_single_source_cgal_program("test_triply_periodic_minimal_surfaces.cpp")
target_link_libraries(test_triply_periodic_minimal_surfaces PUBLIC CGAL::Eigen3_support)
