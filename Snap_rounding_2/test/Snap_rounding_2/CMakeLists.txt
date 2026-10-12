# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.12...3.31)
project(Snap_rounding_2_Tests)

find_package(CGAL REQUIRED)

create_single_source_cgal_program(test_snap_rounding_2.cpp)
create_single_source_cgal_program(test_vertical_slab_snap_rounding.cpp)
create_single_source_cgal_program(test_vertical_slab_snap_rounding_concept.cpp)