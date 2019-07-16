cmake_minimum_required(VERSION 3.1...3.15)
project( Optimal_transportation_reconstruction_2_Tests )

find_package(CGAL QUIET)

if ( CGAL_FOUND )

  include_directories (BEFORE "include")

  create_single_source_cgal_program( "test_basic.cpp" )
  create_single_source_cgal_program( "test_flip_procedure.cpp" )
  create_single_source_cgal_program( "test_output_modules.cpp" )
  create_single_source_cgal_program( "test_quality.cpp" )
  create_single_source_cgal_program( "test_reconstruction_until.cpp" )
  create_single_source_cgal_program( "test_vertex_edge.cpp" )
  create_single_source_cgal_program( "otr2_epeck.cpp" )

else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()

