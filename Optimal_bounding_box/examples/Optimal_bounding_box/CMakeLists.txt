# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.15)
project(Optimal_bounding_box_Examples)

find_package(CGAL REQUIRED)

include(${CGAL_USE_FILE})

find_package(Eigen3 3.1.0 REQUIRED) #(3.1.0 or greater)
include(CGAL_Eigen_support)
if(NOT TARGET CGAL::Eigen_support)
  message(
    STATUS "This project requires the Eigen library, and will not be compiled.")
  return()
endif()

create_single_source_cgal_program("obb_example.cpp")
create_single_source_cgal_program("obb_with_point_maps_example.cpp")
create_single_source_cgal_program("rotated_aabb_tree_example.cpp")

foreach(target obb_example obb_with_point_maps_example
               rotated_aabb_tree_example)
  target_link_libraries(${target} PUBLIC CGAL::Eigen_support)
endforeach()
