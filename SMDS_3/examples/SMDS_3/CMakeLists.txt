# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.1...3.20)

project(SMDS_3_Examples)

# CGAL and its components
find_package(CGAL REQUIRED)

create_single_source_cgal_program("c3t3_example.cpp")
create_single_source_cgal_program("tetrahedron_soup_to_c3t3_example.cpp")
