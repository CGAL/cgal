# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.1...3.15)
project(Periodic_3_mesh_3_Examples)

if(POLICY CMP0074)
  cmake_policy(SET CMP0074 NEW)
endif()

# CGAL and its components
find_package(CGAL REQUIRED COMPONENTS ImageIO)

include(${CGAL_USE_FILE})

# include for local package

# Use Eigen
find_package(Eigen3 3.1.0 REQUIRED) #(3.1.0 or greater)
include(CGAL_Eigen_support)
if(NOT TARGET CGAL::Eigen_support)
  message(
    STATUS "This project requires the Eigen library, and will not be compiled.")
  return()
endif()

# Boost and its components
find_package(Boost)

if(NOT Boost_FOUND)
  message(
    STATUS "This project requires the Boost library, and will not be compiled.")
  return()
endif()

# Creating entries for all C++ files with "main" routine
# ##########################################################

create_single_source_cgal_program("mesh_implicit_shape.cpp")
create_single_source_cgal_program("mesh_implicit_multi_domain.cpp")
create_single_source_cgal_program("mesh_implicit_shape_with_subdomains.cpp")
create_single_source_cgal_program("mesh_implicit_shape_with_optimizers.cpp")
create_single_source_cgal_program("mesh_implicit_shape_with_features.cpp")

foreach(
  target
  mesh_implicit_shape mesh_implicit_multi_domain
  mesh_implicit_shape_with_subdomains mesh_implicit_shape_with_optimizers
  mesh_implicit_shape_with_features)
  target_link_libraries(${target} PUBLIC CGAL::Eigen_support)
endforeach()
