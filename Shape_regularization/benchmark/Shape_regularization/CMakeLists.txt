# Created by the script cgal_create_CMakeLists.
# This is the CMake script for compiling a set of CGAL applications.

project(Shape_regularization_Benchmarks)

cmake_minimum_required(VERSION 3.12...3.29)

find_package(CGAL REQUIRED COMPONENTS Core)

# Find OSQP library and headers.
find_package(OSQP QUIET)
include(CGAL_OSQP_support)
if(TARGET CGAL::OSQP_support)
  message(STATUS "Found OSQP")

  set(osqp_targets benchmark_contours benchmark_qp_segments)
  foreach(osqp_target ${osqp_targets})
    create_single_source_cgal_program("${osqp_target}.cpp")
    if(TARGET ${osqp_target})
      target_link_libraries(${osqp_target} PUBLIC CGAL::OSQP_support)
    endif()
  endforeach()
else()
  message("NOTICE: OSQP was not found. OSQP benchmarks won't be available.")
endif()
