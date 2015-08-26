# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( examples_ )

cmake_minimum_required(VERSION 2.8.10)

find_package(CGAL QUIET COMPONENTS Core )

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )

  include( CGAL_CreateSingleSourceCGALProgram )

  find_package(Eigen3 3.1.0) #(requires 3.1.0 or greater)
  if (NOT EIGEN3_FOUND)
    message(STATUS "This program requires the Eigen library, and will not be compiled.")
  else()
    include( ${EIGEN3_USE_FILE} )
  endif()

  
  include_directories (BEFORE "../include")

  create_single_source_cgal_program( "singular_value_decomposition.cpp" )
  create_single_source_cgal_program( "diagonalize_matrix.cpp" )
  create_single_source_cgal_program( "sparse_solvers.cpp" )

else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()

