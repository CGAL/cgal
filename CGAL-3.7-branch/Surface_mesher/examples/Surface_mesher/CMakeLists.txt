# This is the CMake script for compiling a CGAL application.


project( Surface_mesher_example )

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)
 
if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )  
endif()

include_directories(../../include)

find_package(CGAL QUIET COMPONENTS ImageIO)

if ( CGAL_FOUND AND CGAL_ImageIO_FOUND )

  include( ${CGAL_USE_FILE} )

  include( CGAL_CreateSingleSourceCGALProgram )

  create_single_source_cgal_program( "mesh_a_3d_gray_image.cpp" )
  create_single_source_cgal_program( "mesh_an_implicit_function.cpp" )

else()
  
    message(STATUS "NOTICE: This program requires the CGAL and CGAL ImageIO libraries, and will not be compiled.")
  
endif()

