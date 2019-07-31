# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


cmake_minimum_required(VERSION 3.1...3.15)
project( Periodic_4_hyperbolic_triangulation_2_Benchmarks )


find_package(CGAL QUIET COMPONENTS Core )

if ( CGAL_FOUND )

  	include( ${CGAL_USE_FILE} )

  	include( CGAL_CreateSingleSourceCGALProgram )

    create_single_source_cgal_program( "bench_p4ht2_hyperbolic_vs_euclidean.cpp" ) 
    create_single_source_cgal_program( "bench_p4ht2_insertion.cpp" )
    create_single_source_cgal_program( "bench_p4ht2_remove_dummy_points.cpp" )
    
else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()

