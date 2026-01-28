# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.12...3.31)
project(Optimal_bounding_box_Examples)

find_package(CGAL REQUIRED)

find_package(Eigen3 QUIET)
include(CGAL_Eigen3_support)

if(NOT TARGET CGAL::Eigen3_support)
  message("NOTICE: This project requires the Eigen library, and will not be compiled.")
  return()
endif()

create_single_source_cgal_program("obb_example.cpp")
create_single_source_cgal_program("obb_with_point_maps_example.cpp")
create_single_source_cgal_program("rotated_aabb_tree_example.cpp")

foreach(target obb_example obb_with_point_maps_example rotated_aabb_tree_example)
  target_link_libraries(${target} PRIVATE CGAL::Eigen3_support)
endforeach()
