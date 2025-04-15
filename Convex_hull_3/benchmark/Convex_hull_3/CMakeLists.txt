# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.22)
project(Periodic_2_triangulation_2_Benchs)

find_package(CGAL REQUIRED)

find_package(Eigen3 3.1.0 REQUIRED) #(3.1.0 or greater)
include(CGAL_Eigen3_support)

include_directories(BEFORE "include")

# create a target per cppfile
file(
  GLOB cppfiles
  RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
foreach(cppfile ${cppfiles})
  create_single_source_cgal_program("${cppfile}")
endforeach()

foreach(target bench_ch3_do_intersect)
  target_link_libraries(${target} PUBLIC CGAL::Eigen3_support)
endforeach()