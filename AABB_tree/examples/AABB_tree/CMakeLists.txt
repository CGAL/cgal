# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project(AABB_examples)

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

include_directories(../../include/)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)
 
if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )  
endif()
 
# Find CGAL 
find_package(CGAL COMPONENTS)
#include( ${CGAL_USE_FILE} )

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )

  include( CGAL_CreateSingleSourceCGALProgram )

  create_single_source_cgal_program("AABB_segment_3_example.cpp")
  create_single_source_cgal_program("AABB_triangle_3_example.cpp")
  create_single_source_cgal_program("AABB_polyhedron_edge_example.cpp")
  create_single_source_cgal_program("AABB_polyhedron_facet_distance_example.cpp")
  create_single_source_cgal_program("AABB_polyhedron_facet_intersection_example.cpp")
  create_single_source_cgal_program("AABB_custom_example.cpp")
  create_single_source_cgal_program("AABB_insertion_example.cpp")
  create_single_source_cgal_program("AABB_custom_indexed_triangle_set_example.cpp")
  create_single_source_cgal_program("AABB_custom_triangle_soup_example.cpp")
else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()
