# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.23)
project(Segment_Delaunay_graph_2_example)

find_package(CGAL REQUIRED COMPONENTS Core)

create_single_source_cgal_program("benchmark.cpp")
create_single_source_cgal_program("benchmark_nox.cpp")
create_single_source_cgal_program("double.cpp")
create_single_source_cgal_program("sdg-bench-traits.cpp")
