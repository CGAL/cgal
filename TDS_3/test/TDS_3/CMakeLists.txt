project( TDS_3_Tests )

cmake_minimum_required(VERSION 3.1)

find_package(CGAL QUIET)

if ( CGAL_FOUND )

  include_directories (BEFORE "./include")

  find_package( TBB QUIET )
  create_single_source_cgal_program( "test_triangulation_tds_3.cpp" )
  if(TBB_FOUND)
    CGAL_target_use_TBB(test_triangulation_tds_3)
  endif()
else()
  message(STATUS "This program requires the CGAL library, and will not be compiled.")
endif()
