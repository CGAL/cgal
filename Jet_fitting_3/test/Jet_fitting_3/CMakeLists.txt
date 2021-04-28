# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.20)
project(Jet_fitting_3_Tests)

find_package(CGAL REQUIRED)

# use Eigen
find_package(Eigen3 3.1.0) #(requires 3.1.0 or greater)
include(CGAL_Eigen3_support)
if(TARGET CGAL::Eigen3_support)
  create_single_source_cgal_program("blind_1pt.cpp")
  target_link_libraries(blind_1pt PUBLIC CGAL::Eigen3_support)
else()
  message(
    STATUS
      "NOTICE: This program requires Eigen 3.1 (or greater) and will not be compiled."
  )
endif()
