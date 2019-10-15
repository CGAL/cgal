# Created by the script cgal_create_cmake_script.
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.15)
project(Shape_detection_Benchmarks)

set(CMAKE_CXX_STANDARD 11)

find_package(CGAL QUIET COMPONENTS Core)
if(CGAL_FOUND)

  include(${CGAL_USE_FILE})
  include(CGAL_CreateSingleSourceCGALProgram)

  # Use Eigen.
  find_package(Eigen3 3.1.0 QUIET) # (3.1.0 or greater)
  if(EIGEN3_FOUND)
    include(${EIGEN3_USE_FILE})

    create_single_source_cgal_program(
      "benchmark_region_growing_on_point_set_2.cpp")
    create_single_source_cgal_program(
      "benchmark_region_growing_on_point_set_3.cpp")
  endif()

else()
  message(WARNING 
  "This program requires the CGAL library, and will not be compiled.")
endif()
