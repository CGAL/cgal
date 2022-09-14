# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.23)
project(Triangulation_benchmark)

find_package(CGAL REQUIRED COMPONENTS Core)

find_package(Eigen3 3.1.0 QUIET)
include(CGAL_Eigen3_support)
if(TARGET CGAL::Eigen3_support)
  include_directories(BEFORE "include")

  create_single_source_cgal_program("delaunay.cpp")
  target_link_libraries(delaunay PUBLIC CGAL::Eigen3_support)
  create_single_source_cgal_program("Td_vs_T2_and_T3.cpp")
  target_link_libraries(Td_vs_T2_and_T3 PUBLIC CGAL::Eigen3_support)
else()
  message("NOTICE: Executables in this directory require Eigen 3.1 (or greater), and will not be compiled.")
endif()
