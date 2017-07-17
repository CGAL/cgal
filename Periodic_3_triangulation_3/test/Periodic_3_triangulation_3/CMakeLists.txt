# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Periodic_3_triangulation_3_Tests )

cmake_minimum_required(VERSION 2.8.10)

find_package(CGAL QUIET COMPONENTS Core )

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )

  include( CGAL_CreateSingleSourceCGALProgram )

  include_directories (BEFORE "../../include")

    include_directories (BEFORE "include")

  create_single_source_cgal_program( "test_p3rt3_as_p3t3.cpp" )
  create_single_source_cgal_program( "test_p3rt3_insert_remove_point_set.cpp" )
  create_single_source_cgal_program( "test_periodic_3_alpha_shape_3.cpp" )
  create_single_source_cgal_program( "test_periodic_3_delaunay_3.cpp" )
  create_single_source_cgal_program( "test_periodic_3_delaunay_hierarchy_3.cpp" )
  create_single_source_cgal_program( "test_periodic_3_offset_3.cpp" )
  create_single_source_cgal_program( "test_periodic_3_regular_triangulation_3.cpp" )
  create_single_source_cgal_program( "test_periodic_3_regular_triangulation_traits_3.cpp" )
  create_single_source_cgal_program( "test_periodic_3_triangulation_3.cpp" )
  create_single_source_cgal_program( "test_periodic_3_triangulation_tds.cpp" )
  create_single_source_cgal_program( "test_periodic_3_triangulation_traits_C_3.cpp" )
  create_single_source_cgal_program( "test_periodic_3_triangulation_traits_H_3.cpp" )
  create_single_source_cgal_program( "test_periodic_3_triangulation_traits_SC_3.cpp" )
  create_single_source_cgal_program( "test_periodic_3_triangulation_traits_SH_3.cpp" )

else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()

