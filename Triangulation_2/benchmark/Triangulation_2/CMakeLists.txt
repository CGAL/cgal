# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Triangulation_2_test )

cmake_minimum_required(VERSION 2.6.2)
cmake_policy(VERSION 2.8.4)

find_package(CGAL QUIET COMPONENTS Core )

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )

  include( CGAL_CreateSingleSourceCGALProgram )

  include_directories (BEFORE ../../include)

  include_directories (BEFORE ../../../../experimental-packages/Triangulation_2-unrecursive/include)

  create_single_source_cgal_program( "Triangulation_benchmark_2.cpp" )

else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()

