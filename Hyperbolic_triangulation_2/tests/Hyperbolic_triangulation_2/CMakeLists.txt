# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Hyperbolic_triangulation_2_tests )

cmake_minimum_required(VERSION 2.8.10)

find_package(CGAL QUIET COMPONENTS Core )

if ( CGAL_FOUND )
	
  	include( ${CGAL_USE_FILE} )

  	include( CGAL_CreateSingleSourceCGALProgram )

  	create_single_source_cgal_program( "ht2_test_clear.cpp" )
  	create_single_source_cgal_program( "ht2_test_remove.cpp" )
  	create_single_source_cgal_program( "ht2_test_swap.cpp" )
  	create_single_source_cgal_program( "ht2_test_copy.cpp" )
  	#create_single_source_cgal_program( "ht2_test_equality.cpp" )

else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()

