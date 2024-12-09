# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.12...3.29)
project(Kernel_23_benchmark)

find_package(CGAL QUIET OPTIONAL_COMPONENTS Core)

create_single_source_cgal_program( "cmp_epeck_points.cpp" )
