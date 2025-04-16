# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.23)
project(Ball_merge_surface_reconstruction_Examples)

find_package(CGAL REQUIRED)

create_single_source_cgal_program("ball_merge_reconstruction.cpp")

find_package(TBB QUIET)
include(CGAL_TBB_support)
if(TARGET CGAL::TBB_support)
target_link_libraries(ball_merge_reconstruction PUBLIC CGAL::TBB_support)
else()
  message(STATUS "NOTICE: Intel TBB was not found. Parallel code will not be used.")
endif()
