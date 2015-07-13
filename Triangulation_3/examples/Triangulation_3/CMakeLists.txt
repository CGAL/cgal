# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Triangulation_3_example )

cmake_minimum_required(VERSION 2.8.11)


find_package(CGAL QUIET COMPONENTS Core )

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )
  
  find_package( TBB QUIET )

  if( TBB_FOUND )
    include(${TBB_USE_FILE})
    list(APPEND CGAL_3RD_PARTY_LIBRARIES ${TBB_LIBRARIES})
  endif()

  include( CGAL_CreateSingleSourceCGALProgram )

  include_directories (BEFORE "../../include")

  create_single_source_cgal_program( "adding_handles_3.cpp" )
  create_single_source_cgal_program( "color.cpp" )
  create_single_source_cgal_program( "fast_location_3.cpp" )
  create_single_source_cgal_program( "find_conflicts_3.cpp" )
  create_single_source_cgal_program( "info_insert_with_pair_iterator.cpp" )
  create_single_source_cgal_program( "info_insert_with_pair_iterator_regular.cpp" )
  create_single_source_cgal_program( "info_insert_with_transform_iterator.cpp" )
  create_single_source_cgal_program( "info_insert_with_zip_iterator.cpp" )
  create_single_source_cgal_program( "linking_2d_and_3d.cpp" )
  create_single_source_cgal_program( "sequential_parallel.cpp" )
  create_single_source_cgal_program( "parallel_insertion_in_delaunay_3.cpp" )
  create_single_source_cgal_program( "parallel_insertion_and_removal_in_regular_3.cpp" )
  create_single_source_cgal_program( "regular_3.cpp" )
  create_single_source_cgal_program( "simple_triangulation_3.cpp" )
  create_single_source_cgal_program( "simplex.cpp" )
  create_single_source_cgal_program( "tds.cpp" )

else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()

