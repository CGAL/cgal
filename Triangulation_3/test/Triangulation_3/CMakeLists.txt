# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Triangulation_3_Tests )

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

  create_single_source_cgal_program( "test_delaunay_3.cpp" )
  create_single_source_cgal_program( "test_delaunay_hierarchy_3.cpp" )
  create_single_source_cgal_program( "test_delaunay_hierarchy_3_old.cpp" )
  create_single_source_cgal_program( "test_regular_3.cpp" )
  create_single_source_cgal_program( "test_regular_as_delaunay_3.cpp" )
  create_single_source_cgal_program( "test_regular_insert_range_with_info.cpp" )
  create_single_source_cgal_program( "test_regular_remove_3.cpp" )
  create_single_source_cgal_program( "test_regular_traits_3.cpp" )
  create_single_source_cgal_program( "test_RT_cell_base_with_weighted_circumcenter_3.cpp" )
  create_single_source_cgal_program( "test_robust_weighted_circumcenter.cpp" )
  create_single_source_cgal_program( "test_simplex_3.cpp" )
  create_single_source_cgal_program( "test_static_filters.cpp" )
  create_single_source_cgal_program( "test_triangulation_3.cpp" )

else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()

