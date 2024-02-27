cmake_minimum_required(VERSION 3.1...3.23)
project(Ball_merge_surface_reconstruction_Tests)

find_package(CGAL REQUIRED)

include_directories( ../../include)
create_single_source_cgal_program("ball_merge_reconstruction.cpp")
create_single_source_cgal_program("ball_merge_reconstruction2.cpp")

find_package(TBB QUIET)
include(CGAL_TBB_support)
if(TARGET CGAL::TBB_support)
target_link_libraries(ball_merge_reconstruction PUBLIC CGAL::TBB_support)
target_link_libraries(ball_merge_reconstruction2 PUBLIC CGAL::TBB_support)
else()
  message(STATUS "NOTICE: Intel TBB was not found. Parallel code will not be used.")
endif()
