# Created by the script cgal_create_cmake_script.
# This is the CMake script for compiling a CGAL application.

project(Shape_detection_Benchmarks)

cmake_minimum_required(VERSION 2.8.10)
set(CMAKE_CXX_STANDARD 11)

find_package(CGAL QUIET COMPONENTS Core)
if(CGAL_FOUND)

  include(${CGAL_USE_FILE})
  include(CGAL_CreateSingleSourceCGALProgram)

  create_single_source_cgal_program(
    "benchmark_region_growing_on_point_set_2.cpp")
  create_single_source_cgal_program(
    "benchmark_region_growing_on_point_set_3.cpp")

else()
  message(WARNING 
  "This program requires the CGAL library, and will not be compiled.")
endif()
