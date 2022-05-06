# Created by the script cgal_create_cmake_script.
# This is the CMake script for compiling a CGAL application.

project(Weights_Examples)

cmake_minimum_required(VERSION 3.1...3.22)

find_package(CGAL REQUIRED COMPONENTS Core)

create_single_source_cgal_program("weights.cpp")
create_single_source_cgal_program("coordinates_one_query.cpp")
create_single_source_cgal_program("coordinates_multiple_queries.cpp")
create_single_source_cgal_program("projection_traits.cpp")
create_single_source_cgal_program("custom_traits.cpp")
create_single_source_cgal_program("convergence.cpp")

find_package(Eigen3 3.1.0 QUIET) # (3.1.0 or greater)
include(CGAL_Eigen3_support)
if(TARGET CGAL::Eigen3_support)
  create_single_source_cgal_program("weighted_laplacian.cpp")
  target_link_libraries(weighted_laplacian PUBLIC CGAL::Eigen3_support)
else()
  message(NOTICE "Several examples require the Eigen library, and will not be compiled.")
endif()
