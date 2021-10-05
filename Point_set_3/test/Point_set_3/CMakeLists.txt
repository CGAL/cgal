# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.1...3.20)
project(Point_set_3_Tests)

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

# Creating entries for all C++ files with "main" routine
# ##########################################################

create_single_source_cgal_program("point_set_test.cpp")
create_single_source_cgal_program("point_set_test_join.cpp")
create_single_source_cgal_program("test_deprecated_io_ps.cpp")

#Use LAS
#disable if MSVC 2017
if(NOT MSVC_VERSION OR (MSVC_VERSION GREATER_EQUAL 1919 AND  MSVC_VERSION LESS 1910))
  find_package(LASLIB)
  include(CGAL_LASLIB_support)
  if (TARGET CGAL::LASLIB_support)
    target_link_libraries(test_deprecated_io_ps PUBLIC CGAL::LASLIB_support)
  else()
    message(STATUS "NOTICE : the LAS reader test requires LASlib and will not be compiled.")
  endif()
else()
  message(STATUS "NOTICE : the LAS reader does not work with Visual Studio 2017.")
endif()
