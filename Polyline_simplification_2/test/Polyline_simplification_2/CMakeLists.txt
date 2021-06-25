# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.1...3.20)
project(Polyline_simplification_2_Tests)

# CGAL and its components
find_package(CGAL REQUIRED)

# Boost and its components
find_package(Boost)

if(NOT Boost_FOUND)

  message(
    STATUS "This project requires the Boost library, and will not be compiled.")

  return()

endif()

# include for local directory

# include for local package

# Creating entries for all .cpp/.C files with "main" routine
# ##########################################################

create_single_source_cgal_program( "issue-5774.cpp" )
create_single_source_cgal_program( "simplify_polygon_test.cpp" )

create_single_source_cgal_program( "simplify_polyline_with_duplicate_points.cpp" )
