# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.12...3.29)
project(Stream_support)

# CGAL and its components
find_package(CGAL REQUIRED)

create_single_source_cgal_program("read_doubles.cpp")
create_single_source_cgal_program("read_points.cpp")
