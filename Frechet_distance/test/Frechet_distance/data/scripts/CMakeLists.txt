# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.14)
project( Frechet_distance_naive )

find_package(CGAL REQUIRED)

create_single_source_cgal_program("Compute_naive_Frechet_distance_3.cpp")

