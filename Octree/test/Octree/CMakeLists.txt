# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.14)
project( Octree )

find_package(CGAL REQUIRED QUIET OPTIONAL_COMPONENTS Core )

create_single_source_cgal_program( "octree_test.cpp" )
create_single_source_cgal_program( "test_octree_equality.cpp" )
create_single_source_cgal_program( "test_node_index.cpp" )

