# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.23)
project(Orthtree_Tests)

find_package(CGAL REQUIRED OPTIONAL_COMPONENTS Core)

create_single_source_cgal_program("test_octree_equality.cpp")
create_single_source_cgal_program("test_octree_refine.cpp")
create_single_source_cgal_program("test_octree_grade.cpp")
create_single_source_cgal_program("test_octree_locate.cpp")
create_single_source_cgal_program("test_octree_bbox.cpp")
create_single_source_cgal_program("test_octree_nearest_neighbor.cpp")
create_single_source_cgal_program("test_octree_traverse.cpp")
create_single_source_cgal_program("test_octree_intersecting.cpp")
create_single_source_cgal_program("test_octree_copy_move_constructors.cpp")
create_single_source_cgal_program("test_octree_kernels.cpp")

create_single_source_cgal_program("test_node_index.cpp")
create_single_source_cgal_program("test_node_adjacent.cpp")
