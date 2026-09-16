# Created by the script cgal_create_cmake_script.
# This is the CMake script for compiling a CGAL application.
cmake_minimum_required(VERSION 3.12...3.31)
project(Barycentric_coordinates_2_Benchmarks)

find_package(CGAL REQUIRED COMPONENTS Core)

create_single_source_cgal_program("benchmark_segment_coordinates.cpp")
create_single_source_cgal_program("benchmark_triangle_coordinates.cpp")
create_single_source_cgal_program("benchmark_polygon_4_vertices.cpp")
create_single_source_cgal_program("benchmark_polygon_16_vertices.cpp")
create_single_source_cgal_program("benchmark_polygon_100_vertices.cpp")
create_single_source_cgal_program("benchmark_mv_34_vertices.cpp")

find_package(Eigen3 QUIET)
include(CGAL_Eigen3_support)

if(TARGET CGAL::Eigen3_support)
  create_single_source_cgal_program("benchmark_hm_4_vertices.cpp")
  target_link_libraries(benchmark_hm_4_vertices PRIVATE CGAL::Eigen3_support)
  create_single_source_cgal_program("benchmark_hm_n_vertices.cpp")
  target_link_libraries(benchmark_hm_n_vertices PRIVATE CGAL::Eigen3_support)
else()
  message(STATUS "NOTICE: Several benchmarks require the Eigen library, and will not be compiled.")
endif()
