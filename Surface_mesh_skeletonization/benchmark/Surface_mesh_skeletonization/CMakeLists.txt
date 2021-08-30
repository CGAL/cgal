# Created by the script cgal_create_cmake_script_with_options
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.1...3.20)
project(Mean_curvature_skeleton)

#SET(CMAKE_BUILD_TYPE "Debug")
#SET(GCC_COVERAGE_COMPILE_FLAGS "-fprofile-arcs -ftest-coverage")
#SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}" )

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
find_package(Eigen3 3.2.0) #(requires 3.2.0 or greater)
include(CGAL_Eigen3_support)

if(NOT TARGET CGAL::Eigen3_support)
  message(STATUS "NOTICE: Eigen 3.2 (or greater) is not found.")
endif()

# Creating entries for all .cpp/.C files with "main" routine
# ##########################################################

create_single_source_cgal_program("solver_benchmark.cpp")
target_link_libraries(solver_benchmark PUBLIC CGAL::Eigen3_support)
create_single_source_cgal_program("mcf_scale_invariance.cpp")
target_link_libraries(mcf_scale_invariance PUBLIC CGAL::Eigen3_support)
