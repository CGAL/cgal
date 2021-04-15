# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.1...3.20)
project(Surface_mesh_approximation_Benchmarks)

# CGAL and its components
find_package(CGAL REQUIRED)

# include helper file
include(${CGAL_USE_FILE})

# Boost and its components
find_package(Boost REQUIRED)

if(NOT Boost_FOUND)

  message(
    STATUS "This project requires the Boost library, and will not be compiled.")

  return()

endif()

# include for local directory

# include for local package
include_directories(BEFORE ../../include)

# Creating entries for all C++ files with "main" routine
# ##########################################################

include(CGAL_CreateSingleSourceCGALProgram)

create_single_source_cgal_program("vsa_autoinit_timing_benchmark.cpp")

create_single_source_cgal_program("vsa_timing_benchmark.cpp")
