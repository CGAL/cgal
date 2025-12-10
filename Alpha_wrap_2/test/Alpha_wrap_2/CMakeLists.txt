# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.12...3.31)
project(Alpha_wrap_2_Tests)

find_package(CGAL REQUIRED)

# create a target per cppfile
create_single_source_cgal_program("test_alpha_wrap_2_points.cpp")
create_single_source_cgal_program("test_alpha_wrap_2_polylines.cpp")
create_single_source_cgal_program("test_alpha_wrap_2_API.cpp")
