# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.12...3.29)
project(Hash_map)

# CGAL and its components
find_package(CGAL REQUIRED)

create_single_source_cgal_program("hm.cpp")
create_single_source_cgal_program("foreach.cpp")
create_single_source_cgal_program("triangulation.cpp")
create_single_source_cgal_program("polyhedron.cpp")
create_single_source_cgal_program("surface_mesh.cpp")
create_single_source_cgal_program("arrangement.cpp")
