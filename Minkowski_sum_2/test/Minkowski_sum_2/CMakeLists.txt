# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Minkowski_sum_2_Tests )

cmake_minimum_required(VERSION 2.8.10)

# Commented out C++11 for now
# list(FIND CMAKE_CXX_COMPILE_FEATURES cxx_generalized_initializers has_cpp11)
# if (has_cpp11 LESS 0)
#   message(STATUS "NOTICE: These examples requires a C++11 compiler and will not be compiled.")
#   return()
# endif()

# # Use C++11 for this directory and its sub-directories.
# set(CMAKE_CXX_STANDARD 11)
# set(CMAKE_CXX_STANDARD_REQUIRED TRUE)

find_package(CGAL QUIET COMPONENTS Core )

if ( CGAL_FOUND )

  # create a target per cppfile
  file(GLOB cppfiles RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
  foreach(cppfile ${cppfiles})
    create_single_source_cgal_program( "${cppfile}" )
  endforeach()

else()

    message(STATUS "This program requires the CGAL library, and will not be compiled.")

endif()
