# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

project( Spatial_searching )


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
            
find_package(Eigen3 3.1.91) #(requires 3.2.0 or greater)
if (EIGEN3_FOUND)
  include( ${EIGEN3_USE_FILE} )
endif()

if (MSVC)
    # Turn off VC++ warning
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4244")
endif()


# include for local directory

# include for local package
include_directories( BEFORE ../../include )


# Creating entries for all .cpp/.C files with "main" routine
# ##########################################################

include( CGAL_CreateSingleSourceCGALProgram )

create_single_source_cgal_program( "circular_query.cpp" )

create_single_source_cgal_program( "distance_browsing.cpp" )

create_single_source_cgal_program( "iso_rectangle_2_query.cpp" )

create_single_source_cgal_program( "nearest_neighbor_searching.cpp" )

create_single_source_cgal_program( "searching_with_circular_query.cpp" )

create_single_source_cgal_program( "searching_with_point_with_info.cpp" )

create_single_source_cgal_program( "searching_with_point_with_info_inplace.cpp" )

create_single_source_cgal_program( "searching_with_point_with_info_pmap.cpp" )

create_single_source_cgal_program( "searching_surface_mesh_vertices.cpp" )

create_single_source_cgal_program( "searching_polyhedron_vertices.cpp" )

create_single_source_cgal_program( "user_defined_point_and_distance.cpp" )

create_single_source_cgal_program( "using_fair_splitting_rule.cpp" )

create_single_source_cgal_program( "weighted_Minkowski_distance.cpp" )



if (EIGEN3_FOUND)

  create_single_source_cgal_program( "fuzzy_range_query.cpp" )

  create_single_source_cgal_program( "general_neighbor_searching.cpp" )

else()

  message(STATUS "fuzzy_range_query.cpp and general_neighbor_searching.cpp")
  message(STATUS "will not be compiled as they use CGAL::Epick_d which requires the Eigen library.")

endif()
