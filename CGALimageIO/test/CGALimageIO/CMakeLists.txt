# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( CGALImageIO_test )

cmake_minimum_required(VERSION 2.6.2)

include_directories(../../include)

find_package(CGAL QUIET COMPONENTS ImageIO )

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )

  include( CGAL_CreateSingleSourceCGALProgram )

create_single_source_cgal_program( "test_trilinear_interpolation.cpp" )

else()
  
    message(STATUS "NOTICE: This program requires the CGAL library, and will not be compiled.")
  
endif()

