# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.20)
project(Alpha_wrap_3_Benchmark)

find_package(CGAL REQUIRED)

include_directories (BEFORE ../../include ./Quality ./Robustness) # AW3 includes
include_directories (BEFORE ../../../CGAL-Patches/include)

# create a target per cppfile
create_single_source_cgal_program("Performance/performance_benchmark.cpp")
create_single_source_cgal_program("Quality/quality_benchmark.cpp")
create_single_source_cgal_program("Robustness/robustness_benchmark.cpp")
