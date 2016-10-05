# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Polynomial_kernel_Tests )

cmake_minimum_required(VERSION 2.8.10)

find_package(CGAL QUIET COMPONENTS Core )

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )

  include( CGAL_CreateSingleSourceCGALProgram )

  include_directories (BEFORE "../../include")

    include_directories (BEFORE "include")

  create_single_source_cgal_program( "core_test.cpp" )
  create_single_source_cgal_program( "exact_solvers.cpp" )
  create_single_source_cgal_program( "numeric_solvers.cpp" )
  create_single_source_cgal_program( "polynomial.cpp" )

else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()

