# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

project( BGL_polyhedron_3 )


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

find_package( OpenMesh QUIET )

if ( OpenMesh_FOUND )
include( UseOpenMesh )
else()
  message(STATUS "Examples that use OpenMesh will not be compiled.")
endif()

# include for local directory

# include for local package
include_directories( BEFORE ../../include )


# Creating entries for all C++ files with "main" routine
# ##########################################################

include( CGAL_CreateSingleSourceCGALProgram )

create_single_source_cgal_program( "distance.cpp" )

create_single_source_cgal_program( "incident_vertices.cpp" )

create_single_source_cgal_program( "kruskal.cpp" )

create_single_source_cgal_program( "kruskal_with_stored_id.cpp" )

create_single_source_cgal_program( "normals.cpp" )

create_single_source_cgal_program( "range.cpp" )

create_single_source_cgal_program( "transform_iterator.cpp" )



if(OpenMesh_FOUND)
  create_single_source_cgal_program( "polyhedron_2_OpenMesh.cpp" )
  target_link_libraries( polyhedron_2_OpenMesh ${OPENMESH_LIBRARIES} )
endif()

