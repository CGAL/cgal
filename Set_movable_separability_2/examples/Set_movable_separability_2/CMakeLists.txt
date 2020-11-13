# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.15)
project(Set_movable_separability_2_Examples)

list(FIND CMAKE_CXX_COMPILE_FEATURES cxx_generalized_initializers has_cpp11)
if(has_cpp11 LESS 0)
  message(
    STATUS
      "NOTICE: These examples requires a C++11 compiler and will not be compiled."
  )
  return()
endif()

# Use C++11 for this directory and its sub-directories.
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)

find_package(CGAL REQUIRED)

create_single_source_cgal_program("top_edges_single_mold_trans_cast.cpp")
create_single_source_cgal_program(
  "is_pullout_direction_single_mold_trans_cast.cpp")
create_single_source_cgal_program(
  "pullout_directions_single_mold_trans_cast.cpp")
