# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.12...3.31)
project(Surface_mesh_decomposition_Examples)

# CGAL and its components
find_package(CGAL REQUIRED)

find_package(TBB QUIET)
include(CGAL_TBB_support)

create_single_source_cgal_program("approximate_convex_decomposition.cpp")

if(TARGET CGAL::TBB_support)
  target_link_libraries(approximate_convex_decomposition PRIVATE CGAL::TBB_support)
else()
  message(STATUS "NOTICE: Intel TBB was not found. Sequential code will be used.")
endif()
