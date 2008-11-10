# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( CGALimageIO_test )

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)
 
if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )  
endif()

include_directories(../../include)

find_package(CGAL QUIET COMPONENTS ImageIO )

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )

  include( CGAL_CreateSingleSourceCGALProgram )

create_single_source_cgal_program( "test_trilinear_interpolation.cpp" )

else()
  
    message(STATUS "NOTICE: This program requires the CGAL library, and will not be compiled.")
  
endif()

