# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.12...3.31)

project(Filtered_kernel_3_Benchmarks)

# CGAL and its components
find_package(CGAL REQUIRED)

create_single_source_cgal_program("side_of_oriented_sphere_3.cpp")
target_compile_options(side_of_oriented_sphere_3 PRIVATE /arch:AVX2)
