# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.20)
project(Hyperbolic_triangulation_2_benchmark)

find_package(CGAL REQUIRED COMPONENTS Core)

include(${CGAL_USE_FILE})

include(CGAL_CreateSingleSourceCGALProgram)

create_single_source_cgal_program("bench_insertion_with_different_kernels.cpp")
