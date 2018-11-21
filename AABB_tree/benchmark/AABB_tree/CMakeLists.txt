# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

project( AABB_traits_benchmark)


cmake_minimum_required(VERSION 3.1)
  cmake_policy(VERSION 3.1)


if ( COMMAND cmake_policy )

  cmake_policy( SET CMP0003 NEW )  

endif()

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


