# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.12...3.31)
project(Stream_support)

# CGAL and its components
find_package(CGAL REQUIRED)

find_package(FastFloat QUIET)
include(CGAL_FastFloat_support)

# The usage of FASTFLOAT can be overridden by the user,
# e.g. to use the standard library instead of FastFloat.
# add_definitions(-DCGAL_USE_STD_FROM_CHARS)


create_single_source_cgal_program("read_doubles.cpp")
create_single_source_cgal_program("polygon_mesh_io_benchmark.cpp")

if(TARGET CGAL::FastFloat_support)
  message(STATUS "Found FastFloat")
  target_link_libraries(read_doubles PRIVATE CGAL::FastFloat_support)
  target_link_libraries(polygon_mesh_io_benchmark PRIVATE CGAL::FastFloat_support)
endif()
