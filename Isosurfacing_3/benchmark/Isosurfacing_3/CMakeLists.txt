# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.23)
project( Isosurfacing_3_benchmark )

find_package(CGAL REQUIRED)

create_single_source_cgal_program( "benchmark.cpp" )

find_package(TBB)
include(CGAL_TBB_support)
if(TARGET CGAL::TBB_support)
  target_link_libraries(benchmark PUBLIC CGAL::TBB_support)
endif()

target_compile_definitions(benchmark PUBLIC ${SCENARIO})
target_compile_definitions(benchmark PUBLIC ${KERNEL})
target_compile_definitions(benchmark PUBLIC ${ALGO})
target_compile_definitions(benchmark PUBLIC ${TAG})
