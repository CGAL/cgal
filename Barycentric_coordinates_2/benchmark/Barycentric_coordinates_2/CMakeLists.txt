# Created by the script cgal_create_cmake_script.
# This is the CMake script for compiling a CGAL application.

project(Barycentric_coordinates_2_Benchmarks)

cmake_minimum_required(VERSION 2.8.10)
set(CMAKE_CXX_STANDARD 14)

find_package(CGAL QUIET COMPONENTS Core)
if(CGAL_FOUND)

  include(${CGAL_USE_FILE})
  include(CGAL_CreateSingleSourceCGALProgram)

  create_single_source_cgal_program(
    "benchmark_segment_coordinates.cpp")
  create_single_source_cgal_program(
    "benchmark_triangle_coordinates.cpp")
  create_single_source_cgal_program(
    "benchmark_polygon_4_vertices.cpp")
  create_single_source_cgal_program(
    "benchmark_polygon_16_vertices.cpp")
  create_single_source_cgal_program(
    "benchmark_polygon_100_vertices.cpp")
  create_single_source_cgal_program(
    "benchmark_mv_34_vertices.cpp")

  find_package(Eigen3 REQUIRED)
  if(EIGEN3_FOUND)
    include(${EIGEN3_USE_FILE})

    create_single_source_cgal_program(
      "benchmark_hm_4_vertices.cpp")
    create_single_source_cgal_program(
      "benchmark_hm_n_vertices.cpp")

  else()
    message(WARNING
    "Harmonic coordinates require the Eigen library, and will not be compiled.")
  endif()
else()
  message(WARNING
  "This program requires the CGAL library, and will not be compiled.")
endif()
