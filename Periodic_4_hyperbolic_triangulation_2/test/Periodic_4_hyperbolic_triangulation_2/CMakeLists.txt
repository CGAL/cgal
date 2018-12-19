# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Periodic_4_hyperbolic_triangulation_2_Tests )

cmake_minimum_required(VERSION 2.8.10)

find_package(CGAL QUIET COMPONENTS Core )

if ( CGAL_FOUND )

  	include( ${CGAL_USE_FILE} )

  	include( CGAL_CreateSingleSourceCGALProgram )

  	include_directories (BEFORE "../../include" $HOME/include)

  	create_single_source_cgal_program( "test_p4ht2_dummy_points.cpp" )
  	create_single_source_cgal_program( "test_p4ht2_locate.cpp" )
  	create_single_source_cgal_program( "test_p4ht2_insertion.cpp" )
    create_single_source_cgal_program( "test_p4ht2_intersections.cpp" )
    create_single_source_cgal_program( "test_p4ht2_insertion_euclidean.cpp" )
  	create_single_source_cgal_program( "test_p4ht2_removal.cpp" )
    create_single_source_cgal_program( "test_p4ht2_remove_dummy_points.cpp" )
    create_single_source_cgal_program( "test_p4ht2_complex.cpp" )
  	create_single_source_cgal_program( "test_p4ht2_matrix.cpp" )
    create_single_source_cgal_program( "test_p4ht2_words.cpp" )
    create_single_source_cgal_program( "test_p4ht2_triangulation.cpp" )
    create_single_source_cgal_program( "test_p4ht2_hyperbolic_vs_euclidean.cpp" )
    create_single_source_cgal_program( "test_p4ht2_construct_point_2.cpp" )
	

else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()

