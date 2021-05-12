cmake_minimum_required(VERSION 3.1...3.20)
project(Compact_container_benchmark)

find_package(CGAL REQUIRED)

find_package(TBB)
include(CGAL_TBB_support)
create_single_source_cgal_program("cc_benchmark.cpp")
if(TARGET CGAL::TBB_support)
  target_link_libraries(cc_benchmark PUBLIC CGAL::TBB_support)
endif()
