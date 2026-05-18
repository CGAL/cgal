cmake_minimum_required(VERSION 3.12...3.31)

project(Weights_Examples)

find_package(CGAL REQUIRED COMPONENTS Core)

create_single_source_cgal_program("weights.cpp")
create_single_source_cgal_program("coordinates_one_query.cpp")
create_single_source_cgal_program("coordinates_multiple_queries.cpp")
create_single_source_cgal_program("projection_traits.cpp")
create_single_source_cgal_program("custom_traits.cpp")
create_single_source_cgal_program("convergence.cpp")

find_package(Eigen3 QUIET)
include(CGAL_Eigen3_support)

if(TARGET CGAL::Eigen3_support)
  create_single_source_cgal_program("weighted_laplacian.cpp")
  target_link_libraries(weighted_laplacian PRIVATE CGAL::Eigen3_support)
else()
  message(STATUS "NOTICE: The example 'weighted_laplacian' requires the Eigen library, and will not be compiled.")
endif()
