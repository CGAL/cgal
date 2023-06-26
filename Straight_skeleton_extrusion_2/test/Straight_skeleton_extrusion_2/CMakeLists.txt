# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.23)
project(Straight_skeleton_extrusion_2_Tests)

find_package(CGAL REQUIRED COMPONENTS Qt6 Core)

include_directories(BEFORE "include")

# create a target per cppfile
file(
  GLOB cppfiles
  RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
foreach(cppfile ${cppfiles})
  create_single_source_cgal_program("${cppfile}")
endforeach()

if(CGAL_Qt6_FOUND)
  target_link_libraries(test_sls_extrude PUBLIC CGAL::CGAL_Basic_viewer)
endif()
