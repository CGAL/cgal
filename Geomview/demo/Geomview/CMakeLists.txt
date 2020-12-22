# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.15)
project(Geomview_Demo)

if(NOT POLICY CMP0070 AND POLICY CMP0053)
  # Only set CMP0053 to OLD with CMake<3.10, otherwise there is a warning.
  cmake_policy(SET CMP0053 OLD)
endif()

set_property(DIRECTORY PROPERTY CGAL_NO_TESTING TRUE)

find_package(CGAL REQUIRED)

create_single_source_cgal_program("input.cpp")
create_single_source_cgal_program("kernel.cpp")
create_single_source_cgal_program("gv_terrain.cpp")
