# Created by the script cgal_create_cmake_script_with_options
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.1...3.23)
project(Triangulation_apps)

# CGAL and its components
find_package(CGAL REQUIRED)

find_package(Eigen3 3.1.0 QUIET)
include(CGAL_Eigen3_support)
if(NOT TARGET CGAL::Eigen3_support)
  message("NOTICE: Applications require Eigen 3.1 (or greater), and will not be compiled")
  return()
endif()

# include for local directory
include_directories(BEFORE include)

create_single_source_cgal_program("points_to_RT_to_off.cpp")
target_link_libraries(points_to_RT_to_off PUBLIC CGAL::Eigen3_support)
create_single_source_cgal_program("points_to_DT_to_off.cpp")
target_link_libraries(points_to_DT_to_off PUBLIC CGAL::Eigen3_support)
