# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.12...3.31)
project(Straight_skeleton_3_Tests)

include_directories(BEFORE include)

find_package(CGAL REQUIRED)

find_package(ortools CONFIG REQUIRED)

# create a target per cppfile
file(
  GLOB cppfiles
  RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
foreach(cppfile ${cppfiles})
  get_filename_component(target ${cppfile} NAME_WE)
  create_single_source_cgal_program("${cppfile}")
  target_link_libraries("${target}" PRIVATE ortools::ortools)
endforeach()
