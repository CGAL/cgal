# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project( Periodic_3_mesh_3_Tests )

cmake_minimum_required( VERSION 2.8.11 )

find_package( CGAL QUIET )

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )

  include_directories( BEFORE ../../include )

  include( CGAL_CreateSingleSourceCGALProgram )
  create_single_source_cgal_program( "test_implicit_shapes_bunch.cpp" )
  create_single_source_cgal_program( "test_implicit_shapes_with_protection.cpp" )
  create_single_source_cgal_program( "test_triply_periodic_minimal_surfaces.cpp" )

else()

  message(STATUS "This program requires the CGAL library, and will not be compiled.")

endif()

