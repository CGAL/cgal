# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.14)
project( Isosurfacing_3_examples )

find_package(CGAL REQUIRED)

create_single_source_cgal_program( "simple_marching_cube.cpp" )

find_package(TBB)
include(CGAL_TBB_support)
if(TARGET CGAL::TBB_support)
  target_link_libraries(simple_marching_cube PUBLIC CGAL::TBB_support)
endif()
