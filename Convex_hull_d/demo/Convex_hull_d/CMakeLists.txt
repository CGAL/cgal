# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Convex_hull_d_example )

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)
 
if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )  
endif()
 
find_package(CGAL COMPONENTS Core )

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )

  include( CGAL_CreateSingleSourceCGALProgram )

  create_single_source_cgal_program( "chull_3d-geomview-demo.cpp" )
  create_single_source_cgal_program( "chull_dd-geomview-demo.cpp" )
  create_single_source_cgal_program( "chull_dd-runtime.cpp" )
  create_single_source_cgal_program( "delaunay_dd-runtime.cpp" )

else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()

