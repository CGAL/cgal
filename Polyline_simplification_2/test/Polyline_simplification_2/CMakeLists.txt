# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.12...3.29)
project(Polyline_simplification_2_Tests)

# CGAL and its components
find_package(CGAL REQUIRED)

create_single_source_cgal_program( "issue-5774.cpp" )
create_single_source_cgal_program( "simplify_polygon_test.cpp" )
create_single_source_cgal_program( "simplify_polyline_with_duplicate_points.cpp" )
