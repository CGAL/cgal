# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.12...3.29)

project(Triangulation_3)

# CGAL and its components
find_package(CGAL REQUIRED)

# Boost and its components
find_package(Boost REQUIRED)

if(NOT Boost_FOUND)

  message(
    STATUS "This project requires the Boost library, and will not be compiled.")

  return()

endif()

# include for local directory

# include for local package

# Creating entries for all C++ files with "main" routine
# ##########################################################

create_single_source_cgal_program("ocean.cpp")
create_single_source_cgal_program("incident_edges.cpp")
create_single_source_cgal_program("simple_2.cpp")
create_single_source_cgal_program("simple.cpp")
create_single_source_cgal_program("Triangulation_benchmark_3.cpp")
create_single_source_cgal_program("segment_traverser_benchmark.cpp" )

find_package(benchmark QUIET)
if(NOT TARGET benchmark::benchmark)
  message(STATUS "NOTICE: Some benchmarks require the Google benchmark library, and will not be compiled.")
  return()
endif()

find_package(TBB REQUIRED)
include(CGAL_TBB_support)
if(TARGET CGAL::TBB_support)
  create_single_source_cgal_program("DT3_benchmark_with_TBB.cpp")
  target_link_libraries(DT3_benchmark_with_TBB PRIVATE benchmark::benchmark CGAL::TBB_support)

  add_executable(DT3_benchmark_with_TBB_CCC_approximate_size DT3_benchmark_with_TBB.cpp)
  target_compile_definitions(DT3_benchmark_with_TBB_CCC_approximate_size
                             PRIVATE CGAL_CONCURRENT_COMPACT_CONTAINER_APPROXIMATE_SIZE)
  target_link_libraries(DT3_benchmark_with_TBB_CCC_approximate_size
                        PRIVATE CGAL::CGAL benchmark::benchmark CGAL::TBB_support)
else()
  message(STATUS "NOTICE: Some benchmarks require the TBB library, and will not be compiled.")
endif()
