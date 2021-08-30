cmake_minimum_required(VERSION 3.1...3.20)
project(TDS_3_Tests)

find_package(CGAL REQUIRED)

include_directories(BEFORE "./include")

find_package(TBB QUIET)
include(CGAL_TBB_support)
create_single_source_cgal_program("test_triangulation_tds_3.cpp")
create_single_source_cgal_program("test_io_tds3.cpp")
if(TARGET CGAL::TBB_support)
  target_link_libraries(test_triangulation_tds_3 PUBLIC CGAL::TBB_support)
endif()
