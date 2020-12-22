# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.1...3.14)

project(Tetrahedral_remeshing_Examples)

# CGAL and its components
find_package(CGAL REQUIRED)

# Boost and its components
find_package(Boost REQUIRED)
if(NOT Boost_FOUND)
  message(
    STATUS "This project requires the Boost library, and will not be compiled.")
  return()
endif()

# Use Eigen for Mesh_3
find_package(Eigen3 3.1.0 REQUIRED) #(3.1.0 or greater)
include(CGAL_Eigen_support)
find_package(TBB QUIET)
include(CGAL_TBB_support)

# Concurrent Mesh_3
option(CGAL_ACTIVATE_CONCURRENT_MESH_3 "Activate parallelism in Mesh_3" OFF)
if(CGAL_ACTIVATE_CONCURRENT_MESH_3 OR ENV{CGAL_ACTIVATE_CONCURRENT_MESH_3})
  add_definitions(-DCGAL_CONCURRENT_MESH_3)
  find_package(TBB REQUIRED)
  include(CGAL_TBB_support)
endif()

# Creating entries for all C++ files with "main" routine
# ##########################################################
create_single_source_cgal_program( "tetrahedral_remeshing_example.cpp" )
create_single_source_cgal_program( "tetrahedral_remeshing_with_features.cpp")
create_single_source_cgal_program( "tetrahedral_remeshing_of_one_subdomain.cpp")
create_single_source_cgal_program( "tetrahedral_remeshing_from_mesh.cpp")
create_single_source_cgal_program( "mesh_and_remesh_polyhedral_domain_with_features.cpp" )
target_link_libraries(mesh_and_remesh_polyhedral_domain_with_features PUBLIC CGAL::Eigen_support)
if(CGAL_ACTIVATE_CONCURRENT_MESH_3 AND TARGET CGAL::TBB_support)
  target_link_libraries(mesh_and_remesh_polyhedral_domain_with_features PRIVATE CGAL::TBB_support)
endif()

