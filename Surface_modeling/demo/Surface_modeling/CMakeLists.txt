# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Surface_modeling_demo )

cmake_minimum_required(VERSION 2.8.11)
if(POLICY CMP0043)
  cmake_policy(SET CMP0043 OLD)
endif()

find_package(CGAL QUIET COMPONENTS Core )

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )

  find_package(Eigen3 3.1.91) #(requires 3.2.0 or greater)
  if (EIGEN3_FOUND)
    include( ${EIGEN3_USE_FILE} )
    include( CGAL_CreateSingleSourceCGALProgram )

    include_directories (BEFORE "../../include")
    create_single_source_cgal_program( "deform_mesh_for_botsch08_format.cpp" )

  else()
    message(STATUS "NOTICE: This program requires the Eigen library, version 3.2 or later and will not be compiled.")
  endif()
else()
  
    message(STATUS "NOTICE: This program requires the CGAL library, and will not be compiled.")
  
endif()

