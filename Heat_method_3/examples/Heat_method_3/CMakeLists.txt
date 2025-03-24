# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.12...3.31)
project(Heat_method_3_Examples)

# CGAL and its components
find_package(CGAL REQUIRED)

find_package(Eigen3 3.3.0 QUIET)
include(CGAL_Eigen3_support)
if(NOT TARGET CGAL::Eigen3_support)
  message("NOTICE: These examples require the Eigen library (3.3 or greater), and will not be compiled.")
  return()
endif()

# include for local directory
include_directories(BEFORE include)

create_single_source_cgal_program("heat_method.cpp")
target_link_libraries(heat_method PRIVATE CGAL::Eigen3_support)
create_single_source_cgal_program("heat_method_polyhedron.cpp")
target_link_libraries(heat_method_polyhedron PRIVATE CGAL::Eigen3_support)
create_single_source_cgal_program("heat_method_surface_mesh.cpp")
target_link_libraries(heat_method_surface_mesh PRIVATE CGAL::Eigen3_support)
create_single_source_cgal_program("heat_method_surface_mesh_direct.cpp")
target_link_libraries(heat_method_surface_mesh_direct PRIVATE CGAL::Eigen3_support)
