# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.12...3.29)
project(Filtered_kernel_test)

add_executable(bench_simple_comparisons bench_simple_comparisons.cpp)

find_package(CGAL REQUIRED COMPONENTS Core)

add_executable(bench_orientation_3 "orientation_3.cpp")
target_link_libraries(bench_orientation_3 CGAL::CGAL_Core)

add_executable(bench_comparisons "orientation_3.cpp")
target_link_libraries(bench_comparisons CGAL::CGAL_Core)
set_property(
  TARGET bench_comparisons
  APPEND
  PROPERTY COMPILE_DEFINITIONS ONLY_TEST_COMPARISONS)

get_property(
  DEF
  TARGET bench_comparisons
  PROPERTY COMPILE_DEFINITIONS)
