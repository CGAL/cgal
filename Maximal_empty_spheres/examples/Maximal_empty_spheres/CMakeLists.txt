# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.12...3.31)
project(Maximal_empty_spheres_Examples)

# CGAL and its components
find_package(CGAL REQUIRED)
message("CGAL DIR:${CGAL_DIR}")

find_package(Eigen3 3.2.0 QUIET) #(requires 3.2.0 or greater)
include(CGAL_Eigen3_support)
if(TARGET CGAL::Eigen3_support)
  create_single_source_cgal_program("maximal_empty_spheres_2.cpp")
  target_link_libraries(maximal_empty_spheres_2 PUBLIC CGAL::Eigen3_support)

  create_single_source_cgal_program("maximal_empty_spheres_3.cpp")
  target_link_libraries(maximal_empty_spheres_3 PUBLIC CGAL::Eigen3_support)
else()
  message(STATUS "NOTICE: Examples that use Eigen will not be compiled.")
endif()
