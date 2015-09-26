# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( examples_ )

cmake_minimum_required(VERSION 2.8.10)

find_package(CGAL QUIET COMPONENTS Core )

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )

  include( CGAL_CreateSingleSourceCGALProgram )

  include_directories (BEFORE "../include")
  
  # Use Eigen
  find_package(Eigen3 3.1.0) #(requires 3.1.0 or greater)

  if (EIGEN3_FOUND)
    include( ${EIGEN3_USE_FILE} )
    create_single_source_cgal_program( "singular_value_decomposition.cpp" )
    create_single_source_cgal_program( "sparse_solvers.cpp" )
  endif()

  create_single_source_cgal_program( "diagonalize_matrix.cpp" )

else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()

