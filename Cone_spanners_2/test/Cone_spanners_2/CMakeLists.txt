# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.12...3.29)
project(Cone_spanners_2_Tests)

find_package(CGAL REQUIRED COMPONENTS Core)
find_package(LEDA QUIET)

create_single_source_cgal_program("cones_inexact.cpp")
create_single_source_cgal_program("theta_inexact.cpp")
create_single_source_cgal_program("yao_inexact.cpp")

if(CGAL_Core_FOUND OR LEDA_FOUND)
  create_single_source_cgal_program("cones_exact.cpp")
  create_single_source_cgal_program("theta_exact.cpp")
  create_single_source_cgal_program("yao_exact.cpp")
else()
  message("NOTICE: Some tests require the CGAL_Core library (or LEDA), and will not be compiled.")
endif()
