# Created by the script cgal_create_cmake_script_with_options
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.12...3.29)
project(Mean_curvature_skeleton)

# CGAL and its components
find_package(CGAL REQUIRED)

find_package(Eigen3 3.2.0 QUIET) #(requires 3.2.0 or greater)
include(CGAL_Eigen3_support)
if(TARGET CGAL::Eigen3_support)
  create_single_source_cgal_program("solver_benchmark.cpp")
  target_link_libraries(solver_benchmark PUBLIC CGAL::Eigen3_support)
  create_single_source_cgal_program("mcf_scale_invariance.cpp")
  target_link_libraries(mcf_scale_invariance PUBLIC CGAL::Eigen3_support)
else()
  message("NOTICE: This project requires Eigen 3.2.0 (or greater), and will not be compiled.")
endif()
