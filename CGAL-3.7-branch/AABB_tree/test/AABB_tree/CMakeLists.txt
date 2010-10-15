# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project(AABB_tests)

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

include_directories(../../include)


set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)
 
if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )  
endif()
 
# Find CGAL 
find_package(CGAL COMPONENTS)

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )

  include( CGAL_CreateSingleSourceCGALProgram )

  create_single_source_cgal_program("aabb_intersection_triangle_test.cpp")
  create_single_source_cgal_program("aabb_distance_triangle_test.cpp")
  create_single_source_cgal_program("aabb_distance_edge_test.cpp")
  create_single_source_cgal_program("aabb_distance_triangle_hint_test.cpp")
  create_single_source_cgal_program("aabb_naive_vs_tree_triangle_test.cpp")
  create_single_source_cgal_program("aabb_naive_vs_tree_distance_triangle_test.cpp")
  create_single_source_cgal_program("aabb_naive_vs_tree_distance_segment_test.cpp")
  create_single_source_cgal_program("aabb_correctness_triangle_test.cpp")

  create_single_source_cgal_program("bbox_other_do_intersect_test.cpp")
  create_single_source_cgal_program("triangle_other_intersection_test.cpp")

else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()


