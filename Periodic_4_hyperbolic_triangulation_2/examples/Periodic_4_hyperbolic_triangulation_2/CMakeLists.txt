# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Periodic_4_hyperbolic_triangulation_2_Examples )

cmake_minimum_required(VERSION 3.1...3.13)

find_package(CGAL QUIET COMPONENTS Core )

if ( CGAL_FOUND AND CGAL_Core_FOUND)

  	include( ${CGAL_USE_FILE} )

  	include( CGAL_CreateSingleSourceCGALProgram )

  	create_single_source_cgal_program( "p4ht2_example_insertion.cpp" )

else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()

