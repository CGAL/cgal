# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Surface_mesh_ )

cmake_minimum_required(VERSION 2.8.11)


find_package(CGAL QUIET COMPONENTS Core )

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )

  include( CGAL_CreateSingleSourceCGALProgram )

  include_directories (BEFORE "../../include")

  create_single_source_cgal_program( "sm_join.cpp" )
  create_single_source_cgal_program( "sm_aabbtree.cpp" )
  create_single_source_cgal_program( "sm_bgl.cpp" )
  create_single_source_cgal_program( "sm_circulators.cpp" )
  create_single_source_cgal_program( "sm_do_intersect.cpp" )
  create_single_source_cgal_program( "sm_iterators.cpp" )
  create_single_source_cgal_program( "sm_kruskal.cpp" )
  create_single_source_cgal_program( "sm_memory.cpp" )
  create_single_source_cgal_program( "sm_properties.cpp" )

else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()

