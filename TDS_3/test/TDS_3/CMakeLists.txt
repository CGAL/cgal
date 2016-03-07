# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( TDS_3_test )

cmake_minimum_required(VERSION 2.8.11)


find_package(CGAL QUIET)

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )

  find_package( TBB QUIET )

  if( TBB_FOUND )
    include(${TBB_USE_FILE})
    list(APPEND CGAL_3RD_PARTY_LIBRARIES ${TBB_LIBRARIES})
  endif()
  
  include( CGAL_CreateSingleSourceCGALProgram )

  include_directories (BEFORE "../../include")

    include_directories (BEFORE "include")

  create_single_source_cgal_program( "test_triangulation_tds_3.cpp" )

else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()

