# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.12...3.31)
project(Alpha_wrap_2_Examples)

find_package(CGAL REQUIRED)

create_single_source_cgal_program("AABB_2.cpp")

# create a target per cppfile
create_single_source_cgal_program("triangle_soup_wrap_2.cpp")
create_single_source_cgal_program("point_set_wrap_2.cpp")
create_single_source_cgal_program("mixed_inputs_wrap_2.cpp")
