# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.20)
project(AABB_traits_benchmark)

find_package(CGAL REQUIRED OPTIONAL_COMPONENTS Core)

# google benchmark
find_package(benchmark)

if (benchmark_FOUND)
  create_single_source_cgal_program("tree_creation.cpp")
  target_link_libraries(tree_creation benchmark::benchmark)
endif()
create_single_source_cgal_program("test.cpp")
create_single_source_cgal_program("tree_construction.cpp")
