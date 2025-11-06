# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.12...3.29)
project(Periodic_3_mesh_3_Examples)

# CGAL and its components
find_package(CGAL REQUIRED)

# Use Eigen
find_package(Eigen3 QUIET)
include(CGAL_Eigen3_support)

if(NOT TARGET CGAL::Eigen3_support)
  message("NOTICE: This project requires the Eigen library, and will not be compiled.")
  return()
endif()

create_single_source_cgal_program("mesh_implicit_shape.cpp")
create_single_source_cgal_program("mesh_implicit_multi_domain.cpp")
create_single_source_cgal_program("mesh_implicit_shape_with_subdomains.cpp")
create_single_source_cgal_program("mesh_implicit_shape_with_optimizers.cpp")
create_single_source_cgal_program("mesh_implicit_shape_with_features.cpp")
create_single_source_cgal_program("mesh_periodic_polyhedral_domain.cpp")

foreach(
  target
  mesh_implicit_shape
  mesh_implicit_multi_domain
  mesh_implicit_shape_with_subdomains
  mesh_implicit_shape_with_optimizers
  mesh_implicit_shape_with_features
  mesh_periodic_polyhedral_domain)
  target_link_libraries(${target} PUBLIC CGAL::Eigen3_support)
endforeach()
