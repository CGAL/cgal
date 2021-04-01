# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.14)
project(Orthtree_benchmarks)

# TODO: I shouldn't leave this here
set(CMAKE_CXX_STANDARD 14)

find_package(CGAL REQUIRED QUIET OPTIONAL_COMPONENTS Core)

create_single_source_cgal_program("construction.cpp")
create_single_source_cgal_program("nearest_neighbor.cpp")
