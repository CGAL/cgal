# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Surface_modeling_example )

cmake_minimum_required(VERSION 2.6.2)

find_package(CGAL QUIET COMPONENTS Core )

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )

  include( CGAL_CreateSingleSourceCGALProgram )

  include_directories (BEFORE ../../include)

  create_single_source_cgal_program( "hello.cpp" )
  create_single_source_cgal_program( "k-ring.cpp" )
  create_single_source_cgal_program( "k-ring_BGL.cpp" )
  create_single_source_cgal_program( "example.cpp" )

else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()

