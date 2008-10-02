# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( HalfedgeDS_example )

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)
 
find_package(CGAL QUIET COMPONENTS Core )

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )

  include( CreateSingleSourceCGALProgram )

create_single_source_cgal_program( "hds_prog_color.cpp" )
create_single_source_cgal_program( "hds_prog_compact2.cpp" )
create_single_source_cgal_program( "hds_prog_compact.cpp" )
create_single_source_cgal_program( "hds_prog_default.cpp" )
create_single_source_cgal_program( "hds_prog_edge_iterator.cpp" )
create_single_source_cgal_program( "hds_prog_graph2.cpp" )
create_single_source_cgal_program( "hds_prog_graph.cpp" )
create_single_source_cgal_program( "hds_prog_halfedge_iterator.cpp" )
create_single_source_cgal_program( "hds_prog_vector.cpp" )

else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()

