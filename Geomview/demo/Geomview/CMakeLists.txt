# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Geomview_Demo )

cmake_minimum_required(VERSION 2.8.11)
if(POLICY CMP0043)
  cmake_policy(SET CMP0043 OLD)
endif()

find_package(CGAL QUIET COMPONENTS Core )

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )

  include( CGAL_CreateSingleSourceCGALProgram )

  create_single_source_cgal_program( "input.cpp" )
  create_single_source_cgal_program( "kernel.cpp" )
  create_single_source_cgal_program( "gv_terrain.cpp" )

else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()

