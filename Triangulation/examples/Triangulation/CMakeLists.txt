# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Triangulation_example )

cmake_minimum_required(VERSION 2.8.11)


find_package(CGAL QUIET COMPONENTS Core )

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )

  include( CGAL_CreateSingleSourceCGALProgram )
  
  find_package(Eigen3 3.1.0) 
  if (EIGEN3_FOUND) 
    include( ${EIGEN3_USE_FILE} ) 
    include_directories (BEFORE "../../include")

    create_single_source_cgal_program( "barycentric_subdivision.cpp" )
    create_single_source_cgal_program( "delaunay_triangulation.cpp" )
    create_single_source_cgal_program( "triangulation.cpp" )
    create_single_source_cgal_program( "triangulation_data_structure_dynamic.cpp" )
    create_single_source_cgal_program( "triangulation_data_structure_static.cpp" )

  else()
    message(STATUS "NOTICE: Some of the executables in this directory need Eigen 3.1 (or greater) and will not be compiled.")
  endif()

else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()

