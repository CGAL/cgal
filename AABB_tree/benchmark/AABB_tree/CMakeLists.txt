# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.1...3.15)
project( AABB_traits_benchmark)

# CGAL and its components
find_package( CGAL QUIET)
if ( CGAL_FOUND )
  include( ${CGAL_USE_FILE} )
else ()
  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()  

endif()


# Boost and its components
find_package( Boost REQUIRED )
# include for local directory
if ( NOT Boost_FOUND )
  message(STATUS "This project requires the Boost library, and will not be compiled.")
  return()  
endif()

# include for local package
include_directories( BEFORE ../../include )

add_executable (test_ test.cpp)


