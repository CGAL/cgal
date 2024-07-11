# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.12...3.29)
project(Segment_Delaunay_graph_2_example)

find_package(CGAL REQUIRED COMPONENTS Core)

create_single_source_cgal_program("sdg-creation-time.cpp")
create_single_source_cgal_program("benchmark-gen.cpp")
create_single_source_cgal_program("incirc.cpp")
