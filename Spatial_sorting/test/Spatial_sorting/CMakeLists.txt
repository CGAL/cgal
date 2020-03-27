# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


cmake_minimum_required(VERSION 3.1...3.15)
project( Spatial_sorting_Tests )

find_package(CGAL QUIET)

if ( CGAL_FOUND )

  create_single_source_cgal_program( "test_hilbert.cpp" )
  create_single_source_cgal_program( "test_multiscale.cpp" )

  find_package (TBB QUIET)
  if( TBB_FOUND )
    include( CGAL_target_use_TBB )
    CGAL_target_use_TBB( test_hilbert )
  endif( TBB_FOUND )

else()

    message(STATUS "This program requires the CGAL library, and will not be compiled.")

endif()

