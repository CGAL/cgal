# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Surface_reconstruction_3_example )

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)
 
if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )  
endif()
 
include_directories (BEFORE include)
include_directories (BEFORE ../../include)

find_package(CGAL QUIET COMPONENTS Core )

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )
  include( CGAL_CreateSingleSourceCGALProgram )

  find_package(TAUCS QUIET)

  if(TAUCS_FOUND)

    include( ${TAUCS_USE_FILE} )
    create_single_source_cgal_program( "APSS_reconstruction.cpp" )
    create_single_source_cgal_program( "normal_estimation.cpp" )
    create_single_source_cgal_program( "poisson_reconstruction.cpp" )

  else(TAUCS_FOUND)

    message(STATUS "Notice: This directory needs the TAUCS and will not be compiled.")

  endif(TAUCS_FOUND)

else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()
