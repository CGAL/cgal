# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.12...3.31)
project(Straight_skeleton_3_Examples)

find_package(CGAL REQUIRED)

find_package(ortools CONFIG)

if (NOT TARGET ortools::ortools)
  message(STATUS "NOTICE: The OR-Tools library is not found, and the examples will not be compiled.")
  return()
endif()

create_single_source_cgal_program("create_offset_of_non_triangle_polyhedron.cpp")
target_link_libraries(create_offset_of_non_triangle_polyhedron PRIVATE ortools::ortools)

create_single_source_cgal_program("create_offset_polyhedra.cpp")
target_link_libraries(create_offset_polyhedra PRIVATE ortools::ortools)

create_single_source_cgal_program("create_straight_skeleton_3.cpp")
target_link_libraries(create_straight_skeleton_3 PRIVATE ortools::ortools)

create_single_source_cgal_program("create_offset_and_remesh.cpp")
target_link_libraries(create_offset_and_remesh PRIVATE ortools::ortools)
