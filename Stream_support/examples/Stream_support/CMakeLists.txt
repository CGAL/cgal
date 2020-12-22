# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.1...3.15)
project(Stream_support_Examples)

# CGAL and its components
find_package(CGAL REQUIRED)
include(${CGAL_USE_FILE})

find_package(Boost QUIET)

if(NOT Boost_FOUND)

  message(
    STATUS "This project requires the Boost library, and will not be compiled.")

  return()

endif()

create_single_source_cgal_program("Point_WKT.cpp")
create_single_source_cgal_program("Polygon_WKT.cpp")
create_single_source_cgal_program("Linestring_WKT.cpp")
create_single_source_cgal_program("read_WKT.cpp")

create_single_source_cgal_program("read_xml.cpp")
