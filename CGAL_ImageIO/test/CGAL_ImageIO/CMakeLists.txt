# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( CGAL_ImageIO_Tests )

cmake_minimum_required(VERSION 2.8.11)

find_package(CGAL QUIET COMPONENTS ImageIO )

if ( CGAL_FOUND )

  if (WITH_CGAL_ImageIO) 
    create_single_source_cgal_program( "test_trilinear_interpolation.cpp" )
  else()
    message(STATUS "NOTICE: Some tests require the CGAL_ImageIO library, and will not be compiled.")
  endif()

else()
  
    message(STATUS "NOTICE: This program requires the CGAL library, and will not be compiled.")
  
endif()

