# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.14)
project( AABB_traits_benchmark )

find_package(CGAL REQUIRED QUIET OPTIONAL_COMPONENTS Core )

create_single_source_cgal_program( "test.cpp" )
create_single_source_cgal_program( "tree_construction.cpp" )

