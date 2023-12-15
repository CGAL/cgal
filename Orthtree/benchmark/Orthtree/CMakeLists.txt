# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.23)
project(Orthtree_benchmarks)

find_package(CGAL REQUIRED OPTIONAL_COMPONENTS Core)

create_single_source_cgal_program("construction.cpp")
create_single_source_cgal_program("nearest_neighbor.cpp")
