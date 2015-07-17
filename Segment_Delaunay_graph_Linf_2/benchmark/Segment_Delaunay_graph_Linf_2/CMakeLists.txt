# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Segment_Delaunay_graph_2_example )

cmake_minimum_required(VERSION 2.8.11)


find_package(CGAL QUIET COMPONENTS Core )

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )

  include( CGAL_CreateSingleSourceCGALProgram )

  include_directories (BEFORE ../../include)

  create_single_source_cgal_program( "sdg-creation-time.cpp" )
  create_single_source_cgal_program( "benchmark-gen.cpp" )
  create_single_source_cgal_program( "incirc.cpp" )

else()

    message(STATUS "This program requires the CGAL library, and will not be compiled.")

endif()
