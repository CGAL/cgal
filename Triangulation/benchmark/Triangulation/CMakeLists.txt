# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


cmake_minimum_required(VERSION 3.1...3.15)
project( Triangulation_benchmark )



find_package(CGAL QUIET COMPONENTS Core )

if ( CGAL_FOUND )

  include(${CGAL_USE_FILE})

  find_package(Eigen3 3.1.0) 
  if (EIGEN3_FOUND) 
    include( ${EIGEN3_USE_FILE} ) 
    include_directories (BEFORE "include")
    create_single_source_cgal_program( "delaunay.cpp" )
    create_single_source_cgal_program( "Td_vs_T2_and_T3.cpp" )

  else()
    message(STATUS "NOTICE: Some of the executables in this directory need Eigen 3.1 (or greater) and will not be compiled.")
  endif()

else()
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
endif()

