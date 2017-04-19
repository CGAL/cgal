# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

project( BGL_LCC_Examples )


cmake_minimum_required(VERSION 2.8.11)

# CGAL and its components
find_package( CGAL QUIET COMPONENTS  )

if ( NOT CGAL_FOUND )

  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()  

endif()

# include helper file
include( ${CGAL_USE_FILE} )


# Boost and its components
find_package( Boost REQUIRED )

if ( NOT Boost_FOUND )

  message(STATUS "This project requires the Boost library, and will not be compiled.")

  return()  

endif()

# include for local directory
include_directories( BEFORE ../../include )


# Creating entries for all C++ files with "main" routine
# ##########################################################

include( CGAL_CreateSingleSourceCGALProgram )

create_single_source_cgal_program( "distance_lcc.cpp" )

create_single_source_cgal_program( "incident_vertices_lcc.cpp" )

create_single_source_cgal_program( "kruskal_lcc.cpp" )

create_single_source_cgal_program( "normals_lcc.cpp" )

create_single_source_cgal_program( "range_lcc.cpp" )

create_single_source_cgal_program( "transform_iterator_lcc.cpp" )

create_single_source_cgal_program( "copy_lcc.cpp" )

