# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.14)
project( Isosurfacing_3_tests )

find_package(CGAL REQUIRED)

create_single_source_cgal_program( "test.cpp" )

find_package(TBB)
include(CGAL_TBB_support)
if(TARGET CGAL::TBB_support)
  target_link_libraries(test PUBLIC CGAL::TBB_support)
endif()
