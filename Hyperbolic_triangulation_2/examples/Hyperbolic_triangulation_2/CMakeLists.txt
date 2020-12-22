# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.15)
project(Hyperbolic_triangulation_2_Examples)

find_package(CGAL REQUIRED COMPONENTS Core)

include(${CGAL_USE_FILE})

include(CGAL_CreateSingleSourceCGALProgram)

create_single_source_cgal_program("ht2_example.cpp")
create_single_source_cgal_program("ht2_example_color.cpp")
