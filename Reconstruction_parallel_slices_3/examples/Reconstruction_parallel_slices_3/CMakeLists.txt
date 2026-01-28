# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.15)
project(Reconstruction_parallel_slices_3_Examples)

find_package(CGAL REQUIRED)

create_single_source_cgal_program( "conforming-Delaunay_3.cpp" )
create_single_source_cgal_program( "manifold_check.cpp" )
create_single_source_cgal_program( "check_valid_input.cpp" )
create_single_source_cgal_program( "edge_profile.cpp" )
create_single_source_cgal_program( "repair.cpp" )
create_single_source_cgal_program( "extract.cpp" )
create_single_source_cgal_program( "detect_intersecting_segments.cpp" )

