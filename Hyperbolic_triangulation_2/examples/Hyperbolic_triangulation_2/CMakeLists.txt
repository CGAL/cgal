# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.12...3.29)
project(Hyperbolic_triangulation_2_Examples)

find_package(CGAL REQUIRED COMPONENTS Core)
find_package(LEDA QUIET)

if (CGAL_Core_FOUND OR LEDA_FOUND)
  create_single_source_cgal_program("ht2_example.cpp")
  create_single_source_cgal_program("ht2_example_color.cpp")
else()
  message("NOTICE: Examples require CGAL_Core (or LEDA), and will not be compiled.")
endif()
