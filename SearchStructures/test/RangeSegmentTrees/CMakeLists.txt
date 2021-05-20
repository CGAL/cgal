# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.20)
project(RangeSegmentTrees_Tests)

find_package(CGAL REQUIRED)

create_single_source_cgal_program("test_segment_tree_set_2.cpp")
