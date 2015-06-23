# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

project( Spatial_searching )


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

# include for local directory

# include for local package
include_directories( BEFORE ../../include )


# Creating entries for all .cpp/.C files with "main" routine
# ##########################################################

include( CGAL_CreateSingleSourceCGALProgram )

create_single_source_cgal_program( "Building_kd_tree_with_ddim_points.cpp" )

create_single_source_cgal_program( "Building_kd_tree_with_own_pointtype.cpp" )

create_single_source_cgal_program( "Circular_query.cpp" )

create_single_source_cgal_program( "Compare_methods.cpp" )

create_single_source_cgal_program( "Iso_rectangle_2_query.cpp" )

create_single_source_cgal_program( "iso_rectangle_2_query_2.cpp" )

create_single_source_cgal_program( "K_neighbor_search_manhattan_distance_isobox_point.cpp" )

create_single_source_cgal_program( "K_neighbor_search_with_circle.cpp" )

create_single_source_cgal_program( "Orthogonal_incremental_neighbor_search.cpp" )

create_single_source_cgal_program( "Orthogonal_k_neighbor_search.cpp" )

create_single_source_cgal_program( "Range_searching.cpp" )

create_single_source_cgal_program( "rangeVSincremental.cpp" )

create_single_source_cgal_program( "Splitters.cpp" )


