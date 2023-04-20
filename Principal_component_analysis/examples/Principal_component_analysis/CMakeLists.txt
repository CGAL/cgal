# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.23)
project(Principal_component_analysis_Examples)

find_package(CGAL REQUIRED)

# Use Eigen
find_package(Eigen3 3.1.0 QUIET) #(requires 3.1.0 or greater)
include(CGAL_Eigen3_support)
if(NOT TARGET CGAL::Eigen3_support)
  message("NOTICE: This project requires Eigen 3.1 (or greater), and will not be compiled.")
  return()
endif()

# create a target per cppfile
file(
  GLOB cppfiles
  RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
foreach(cppfile ${cppfiles})
  create_single_source_cgal_program("${cppfile}")
  get_filename_component(target ${cppfile} NAME_WE)
  target_link_libraries(${target} PUBLIC CGAL::Eigen3_support)
endforeach()
