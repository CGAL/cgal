# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Ridges_3_test )

cmake_minimum_required(VERSION 2.6.2)
cmake_policy(VERSION 2.8.4)

find_package(CGAL QUIET COMPONENTS Core )

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )

  include( CGAL_CreateSingleSourceCGALProgram )

  # Link with BLAS and LAPACK (required)
  find_package(LAPACK)

  if(LAPACK_FOUND)

    include( ${LAPACK_USE_FILE} )

    create_single_source_cgal_program( "ridge_test.cpp" )

  else(LAPACK_FOUND)
  
    message(STATUS "NOTICE: This program requires LAPACK, and will not be compiled.")

  endif(LAPACK_FOUND)

else()
  
  message(STATUS "NOTICE: This program requires the CGAL library, and will not be compiled.")
  
endif()

