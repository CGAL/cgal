# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.12...3.31)
project(Straight_skeleton_3_Examples)

find_package(CGAL REQUIRED)

find_package(ortools CONFIG REQUIRED)

create_single_source_cgal_program("or_tools_dummy.cpp")
target_link_libraries(or_tools_dummy PRIVATE ortools::ortools)

create_single_source_cgal_program("or_tools_AI.cpp")
target_link_libraries(or_tools_AI PRIVATE ortools::ortools)

create_single_source_cgal_program("create_offset_of_non_triangle_polyhedron.cpp")
target_link_libraries(or_tools_AI PRIVATE ortools::ortools)

create_single_source_cgal_program("create_offset_polyhedra.cpp")
target_link_libraries(or_tools_AI PRIVATE ortools::ortools)

create_single_source_cgal_program("create_straight_skeleton_3.cpp")
target_link_libraries(create_straight_skeleton_3 PRIVATE ortools::ortools)

find_package(Eigen3 QUIET)
include(CGAL_Eigen3_support)

if(TARGET CGAL::Eigen3_support)
  create_single_source_cgal_program("create_offset_and_remesh.cpp")
  target_link_libraries(create_offset_and_remesh PRIVATE CGAL::Eigen3_support ortools::ortools)
else()
  message(STATUS "NOTICE: Several examples require the Eigen library, and will not be compiled.")
endif()