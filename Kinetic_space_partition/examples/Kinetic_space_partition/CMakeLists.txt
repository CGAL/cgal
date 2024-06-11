# Created by the script cgal_create_CMakeLists.
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.1...3.23)

project(Kinetic_space_partition_Examples)

find_package(CGAL REQUIRED)
include(CGAL_CreateSingleSourceCGALProgram)

find_package(Eigen3 3.1.0 REQUIRED)
if(Eigen3_FOUND)
  message(STATUS "Found Eigen")
  include(CGAL_Eigen_support)

  set(targets kinetic_partition)

  foreach(target ${targets})
    create_single_source_cgal_program("${target}.cpp")
    target_link_libraries(${target} PUBLIC CGAL::Eigen_support)
  endforeach()
else()
  message(ERROR "This program requires the Eigen library, and will not be compiled.")
endif()
