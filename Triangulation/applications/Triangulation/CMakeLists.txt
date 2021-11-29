# Created by the script cgal_create_cmake_script_with_options
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.1...3.20)
project(Triangulation_apps)

# CGAL and its components
find_package(CGAL REQUIRED)

# Boost and its components
find_package(Boost REQUIRED)

if(NOT Boost_FOUND)

  message(
    STATUS "This project requires the Boost library, and will not be compiled.")

  return()

endif()

find_package(Eigen3 3.1.0)
include(CGAL_Eigen3_support)

# include for local directory
include_directories(BEFORE include)

# include for local package

# Creating entries for all .cpp/.C files with "main" routine
# ##########################################################

create_single_source_cgal_program("points_to_RT_to_off.cpp")
target_link_libraries(points_to_RT_to_off PUBLIC CGAL::Eigen3_support)
create_single_source_cgal_program("points_to_DT_to_off.cpp")
target_link_libraries(points_to_DT_to_off PUBLIC CGAL::Eigen3_support)
