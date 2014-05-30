# Created by the script cgal_create_cmake_script_with_options
# This is the CMake script for compiling a set of CGAL applications.

project( Surface_mesh_simplification )


cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" VERSION_GREATER 2.6)
  if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
    cmake_policy(VERSION 2.8.4)
  else()
    cmake_policy(VERSION 2.6)
  endif()
endif()

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
find_package( Boost REQUIRED )

if ( NOT Boost_FOUND )

  message(STATUS "This project requires the Boost library, and will not be compiled.")

  return()  

endif()


# Activate concurrency?
option(ACTIVATE_CONCURRENCY
       "Enable concurrency"
       ON)
 
if( ACTIVATE_CONCURRENCY )
  find_package( TBB REQUIRED )
  if( TBB_FOUND )
    include(${TBB_USE_FILE})
    list(APPEND CGAL_3RD_PARTY_LIBRARIES ${TBB_LIBRARIES}) 
  endif()
endif()





# include for local directory

# include for local package
include_directories( BEFORE ../../include )

include_directories (BEFORE "c:/3rdPartyLibs/OpenMesh/src")

add_definitions(
      -D_USE_MATH_DEFINES -DNOMINMAX
  )

link_directories ( c:/3rdPartyLibs/OpenMesh/VC11/Build/lib )


# Creating entries for all .cpp/.C files with "main" routine
# ##########################################################

include( CGAL_CreateSingleSourceCGALProgram )

create_single_source_cgal_program( "tbb.cpp" )

create_single_source_cgal_program( "edge_collapse_constrain_sharp_edges.cpp" )

create_single_source_cgal_program( "edge_collapse_constrained_border_polyhedron.cpp" )

create_single_source_cgal_program( "edge_collapse_enriched_polyhedron.cpp" )

create_single_source_cgal_program( "edge_collapse_polyhedron.cpp" )


create_single_source_cgal_program( "edge_collapse_OpenMesh.cpp" )

target_link_libraries( edge_collapse_OpenMesh debug OpenMeshCored optimized OpenMeshCore )
