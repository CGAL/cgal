# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.14)
project( classical_Frechet_distance )

find_package(CGAL REQUIRED)

create_single_source_cgal_program("Compute_classical_Frechet_distance_3.cpp")
create_single_source_cgal_program("Compute_classical_Frechet_distance_100.cpp")
