# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.12...3.29)
project(Hyperbolic_triangulation_2_benchmark)

find_package(CGAL REQUIRED COMPONENTS Core)

create_single_source_cgal_program("bench_insertion_with_different_kernels.cpp")
