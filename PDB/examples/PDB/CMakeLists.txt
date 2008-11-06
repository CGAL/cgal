# This is the CMake script for compiling a CGAL application.

project( PDB_example ) 

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )  
endif()

find_package(CGAL QUIET COMPONENTS PDB )
include( ${CGAL_USE_FILE} )

if(CGAL_PDB_FOUND)

  include( CGAL_CreateSingleSourceCGALProgram )
  
  create_single_source_cgal_program( "extracting_geometry.cpp" )

else(CGAL_PDB_FOUND)

  message(STATUS "NOTICE: This demo needs the CGAL PDB library, and will not be compiled.")

endif(CGAL_PDB_FOUND)

