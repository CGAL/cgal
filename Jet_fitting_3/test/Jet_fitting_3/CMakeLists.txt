# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Jet_fitting_3_test )

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

find_package(CGAL QUIET COMPONENTS Core )

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )

  include( CreateSingleSourceCGALProgram )

  # Link with BLAS and LAPACK (required)
  include( CGAL_SetupBLAS )
  include( CGAL_SetupLAPACK )

  create_single_source_cgal_program( "blind_1pt.cpp" )

else()

    message(STATUS "This program requires the CGAL library, and will not be compiled.")

endif()

