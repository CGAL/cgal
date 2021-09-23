# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.20)
project(Spatial_sorting_)

find_package(CGAL REQUIRED COMPONENTS Core)

include(${CGAL_USE_FILE})
create_single_source_cgal_program("simple.cpp")
