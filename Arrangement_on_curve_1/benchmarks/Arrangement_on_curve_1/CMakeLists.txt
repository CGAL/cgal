# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.12...3.31)
project(Arrangement_on_curve_1_Benchmarks)

find_package(CGAL REQUIRED)

find_package(Boost COMPONENTS program_options)
if(NOT Boost_PROGRAM_OPTIONS_FOUND)
  message(STATUS "NOTICE: The benchmark requires Boost Program Options, and will not be compiled.")
  return()
endif()

create_single_source_cgal_program("arr_on_circle.cpp")
target_link_libraries(arr_on_circle PRIVATE Boost::program_options)
