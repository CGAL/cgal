# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


cmake_minimum_required(VERSION 3.1...3.15)
project( Ridges_3_Tests )



find_package(CGAL QUIET)

if ( CGAL_FOUND )

  # use either Eigen
  find_package(Eigen3 3.1.0) #(requires 3.1.0 or greater)
  if (EIGEN3_FOUND)
    include( ${EIGEN3_USE_FILE} )
    create_single_source_cgal_program( "ridge_test.cpp" )

  else()
  
    message(STATUS "NOTICE: This program requires Eigen 3.1 (or greater) and will not be compiled.")

  endif()

else()
  
  message(STATUS "NOTICE: This program requires the CGAL library, and will not be compiled.")
  
endif()

