# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.20)
project(Filtered_kernel_test)

add_executable(bench_simple_comparisons bench_simple_comparisons.cpp)

find_package(CGAL REQUIRED COMPONENTS Core)

include(${CGAL_USE_FILE})

add_executable(bench_orientation_3 "orientation_3.cpp")
target_link_libraries(bench_orientation_3 ${CGAL_LIBRARIES}
                      ${CGAL_3RD_PARTY_LIBRARIES})

add_executable(bench_comparisons "orientation_3.cpp")
target_link_libraries(bench_comparisons ${CGAL_LIBRARIES}
                      ${CGAL_3RD_PARTY_LIBRARIES})
set_property(
  TARGET bench_comparisons
  APPEND
  PROPERTY COMPILE_DEFINITIONS ONLY_TEST_COMPARISONS)

get_property(
  DEF
  TARGET bench_comparisons
  PROPERTY COMPILE_DEFINITIONS)
