# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Triangulation_3_Geomview_demos_Demo )

cmake_minimum_required(VERSION 3.0)
if(POLICY CMP0053)
  cmake_policy(SET CMP0053 OLD)
endif()
if(POLICY CMP0043)
  cmake_policy(SET CMP0043 OLD)
endif()

set_directory_properties(PROPERTIES CGAL_NO_TESTING TRUE)

find_package(CGAL QUIET)

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )

  include( CGAL_CreateSingleSourceCGALProgram )

  create_single_source_cgal_program( "Triangulation_3_demo.cpp" )
  create_single_source_cgal_program( "Triangulation_3_color_demo.cpp" )
  create_single_source_cgal_program( "Triangulation_3_remove_demo.cpp" )
  create_single_source_cgal_program( "Triangulation_3_voronoi_demo.cpp" )

else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()

