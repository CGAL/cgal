cmake_minimum_required(VERSION 3.1...3.15)
project( Compact_container_benchmark )


find_package(CGAL QUIET)

if ( CGAL_FOUND )
  find_package( TBB )
  include(CGAL_TBB_support)
  create_single_source_cgal_program( "cc_benchmark.cpp" )
  if(TARGET CGAL::TBB_support)
    target_link_libraries(cc_benchmark PUBLIC CGAL::TBB_support)
  endif()
else()
  message(STATUS "This program requires the CGAL library, and will not be compiled.")
endif()
