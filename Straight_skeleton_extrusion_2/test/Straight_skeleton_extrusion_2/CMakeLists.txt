# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.23)
project(Straight_skeleton_extrusion_2_Tests)

find_package(CGAL REQUIRED COMPONENTS Qt6 Core)
find_package(LEDA QUIET)

include_directories(BEFORE "include")

if (CGAL_Core_FOUND OR LEDA_FOUND)
  create_single_source_cgal_program("test_sls_extrude.cpp")
  if(CGAL_Qt6_FOUND)
    target_link_libraries(test_sls_extrude PUBLIC CGAL::CGAL_Basic_viewer)
  endif()
else()
  message("NOTICE: Tests require CGAL_Core (or LEDA), and will not be compiled.")
endif()
