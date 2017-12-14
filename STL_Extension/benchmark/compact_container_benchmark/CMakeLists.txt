project( Compact_container_benchmark )

cmake_minimum_required(VERSION 3.1)

find_package(CGAL QUIET)

if ( CGAL_FOUND )
  find_package( TBB )
  create_single_source_cgal_program( "cc_benchmark.cpp" )
  if(TBB_FOUND)
    CGAL_target_use_TBB(cc_benchmark)
  endif()
else()
  message(STATUS "This program requires the CGAL library, and will not be compiled.")
endif()
