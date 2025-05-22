# Created by the script cgal_create_cmake_script.
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.12...3.31)

project(Barycentric_coordinates_2_Examples)

find_package(CGAL REQUIRED COMPONENTS Core)

create_single_source_cgal_program("segment_coordinates.cpp")
create_single_source_cgal_program("triangle_coordinates.cpp")
create_single_source_cgal_program("wachspress_coordinates.cpp")
create_single_source_cgal_program("mean_value_coordinates.cpp")
create_single_source_cgal_program("discrete_harmonic_coordinates.cpp")
create_single_source_cgal_program("terrain_height_modeling.cpp")

# this code is deprecated:
create_single_source_cgal_program("deprecated_coordinates.cpp")

find_package(Eigen3 3.1.0 QUIET) # (3.1.0 or greater)
include(CGAL_Eigen3_support)
if(TARGET CGAL::Eigen3_support)
  create_single_source_cgal_program("affine_coordinates.cpp")
  target_link_libraries(affine_coordinates PRIVATE CGAL::Eigen3_support)
  create_single_source_cgal_program("harmonic_coordinates.cpp")
  target_link_libraries(harmonic_coordinates PRIVATE CGAL::Eigen3_support)
  create_single_source_cgal_program("shape_deformation.cpp")
  target_link_libraries(shape_deformation PRIVATE CGAL::Eigen3_support)
else()
  message(STATUS "NOTICE: Several examples require the Eigen library, and will not be compiled.")
endif()
