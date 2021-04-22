# Created by the script cgal_create_cmake_script.
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.20)
project(Shape_detection_Benchmarks)

find_package(CGAL REQUIRED COMPONENTS Core)

include(${CGAL_USE_FILE})
include(CGAL_CreateSingleSourceCGALProgram)

# Use Eigen.
find_package(Eigen3 3.1.0 QUIET) # (3.1.0 or greater)
include(CGAL_Eigen3_support)
if(TARGET CGAL::Eigen3_support)
  create_single_source_cgal_program(
    "benchmark_region_growing_on_point_set_2.cpp")
  target_link_libraries(benchmark_region_growing_on_point_set_2
                        PUBLIC CGAL::Eigen3_support)
  create_single_source_cgal_program(
    "benchmark_region_growing_on_point_set_3.cpp")
  target_link_libraries(benchmark_region_growing_on_point_set_3
                        PUBLIC CGAL::Eigen3_support)
endif()
