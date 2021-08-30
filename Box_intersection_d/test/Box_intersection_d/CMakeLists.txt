# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.20)
project(Box_intersection_d_Tests)

find_package(CGAL REQUIRED)

find_package(TBB)
include(CGAL_TBB_support)

create_single_source_cgal_program("automated_test.cpp")
create_single_source_cgal_program("benchmark_box_intersection.cpp")
create_single_source_cgal_program("random_set_test.cpp")
create_single_source_cgal_program("test_box_grid.cpp")

if(TARGET CGAL::TBB_support)
  target_link_libraries(test_box_grid PUBLIC CGAL::TBB_support)
else()
  message(
    STATUS "NOTICE: Intel TBB was not found. Sequential code will be used.")
endif()
