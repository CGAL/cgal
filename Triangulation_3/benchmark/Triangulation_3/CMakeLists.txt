
project( Triangulation_3_Benchmark )

cmake_minimum_required(VERSION 2.8.11)

find_package(CGAL QUIET)

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )

#  find_package( TBB QUIET )

  include( CGAL_CreateSingleSourceCGALProgram )

  include_directories (BEFORE "../../include")

  create_single_source_cgal_program( "incident_edges.cpp" )
  create_single_source_cgal_program( "simple.cpp" )
  create_single_source_cgal_program( "simple_2.cpp" )
  create_single_source_cgal_program( "Triangulation_benchmark_3.cpp" )
  create_single_source_cgal_program( "segment_traverser_benchmark.cpp" )

else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()

