# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Jet_fitting_3_test )

cmake_minimum_required(VERSION 2.6.2)

find_package(CGAL QUIET COMPONENTS Core )

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )

  include( CGAL_CreateSingleSourceCGALProgram )

  # Link with BLAS and LAPACK (required)
  find_package(LAPACK)

  if(LAPACK_FOUND)

    include( ${LAPACK_USE_FILE} )

  include_directories (BEFORE ../../include)

    create_single_source_cgal_program( "blind_1pt.cpp" )

  else(LAPACK_FOUND)

    message(STATUS "NOTICE: This program requires LAPACK, and will not be compiled.")

  endif(LAPACK_FOUND)

else()

  message(STATUS "NOTICE: This program requires the CGAL library, and will not be compiled.")

endif()

