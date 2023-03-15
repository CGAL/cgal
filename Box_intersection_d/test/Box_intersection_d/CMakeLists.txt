# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.23)
project(Box_intersection_d_Tests)

find_package(CGAL REQUIRED)

create_single_source_cgal_program("automated_test.cpp")
create_single_source_cgal_program("benchmark_box_intersection.cpp")
create_single_source_cgal_program("random_set_test.cpp")
create_single_source_cgal_program("test_box_grid.cpp")
create_single_source_cgal_program("test_Has_member_report.cpp")

find_package(TBB QUIET)
include(CGAL_TBB_support)
if(TARGET CGAL::TBB_support)
  target_link_libraries(test_box_grid PUBLIC CGAL::TBB_support)
else()
  message(STATUS "NOTICE: Intel TBB was not found. Parallel code will not be used.")
endif()
