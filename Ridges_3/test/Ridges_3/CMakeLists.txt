# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Ridges_3_test )

cmake_minimum_required(VERSION 2.8.11)


find_package(CGAL QUIET)

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )

  include( CGAL_CreateSingleSourceCGALProgram )

  # use either Eigen or BLAS/LAPACK
  find_package(Eigen3 3.1.0) #(requires 3.1.0 or greater)
  if (NOT EIGEN3_FOUND)
    find_package(LAPACK)
    if(LAPACK_FOUND)
      include( ${LAPACK_USE_FILE} )
    endif(LAPACK_FOUND)
  else()
    include( ${EIGEN3_USE_FILE} )
  endif()
  
  
  if(EIGEN3_FOUND OR LAPACK_FOUND)

    create_single_source_cgal_program( "ridge_test.cpp" )

  else(EIGEN3_FOUND OR LAPACK_FOUND)
  
    message(STATUS "NOTICE: This program requires either Eigen 3.1 (or greater) or LAPACK, and will not be compiled.")

  endif(EIGEN3_FOUND OR LAPACK_FOUND)

else()
  
  message(STATUS "NOTICE: This program requires the CGAL library, and will not be compiled.")
  
endif()

