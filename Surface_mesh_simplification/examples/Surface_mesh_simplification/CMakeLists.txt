# Created by the script cgal_create_cmake_script_with_options
# This is the CMake script for compiling a set of CGAL applications.

project( Surface_mesh_simplification )


cmake_minimum_required(VERSION 2.8.11)


set( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true )
 
if ( COMMAND cmake_policy )

  cmake_policy( SET CMP0003 NEW )  

endif()

# CGAL and its components
find_package( CGAL QUIET COMPONENTS  )

if ( NOT CGAL_FOUND )

  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()  

endif()

# include helper file
include( ${CGAL_USE_FILE} )


# Boost and its components
find_package( Boost )

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



# Creating entries for all .cpp/.C files with "main" routine
# ##########################################################

include( CGAL_CreateSingleSourceCGALProgram )


create_single_source_cgal_program( "edge_collapse_constrain_sharp_edges.cpp" )

create_single_source_cgal_program( "edge_collapse_constrained_border_polyhedron.cpp" )

create_single_source_cgal_program( "edge_collapse_enriched_polyhedron.cpp" )

create_single_source_cgal_program( "edge_collapse_polyhedron.cpp" )

create_single_source_cgal_program( "edge_collapse_surface_mesh.cpp" )
create_single_source_cgal_program( "edge_collapse_constrained_border_surface_mesh.cpp" )

if(OpenMesh_FOUND)
create_single_source_cgal_program( "edge_collapse_OpenMesh.cpp" )

target_link_libraries( edge_collapse_OpenMesh ${OPENMESH_LIBRARIES} )
endif()
