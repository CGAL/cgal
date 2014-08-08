# Created by the script cgal_create_cmake_script_with_options
# This is the CMake script for compiling a set of CGAL applications.

project( Polyhedron_shortest_path )

SET(CMAKE_CXX_FLAGS "-Wall")

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

find_package( OpenMesh QUIET )

if ( OpenMesh_FOUND )
include( UseOpenMesh )
else()
  message(STATUS "Examples that use OpenMesh will not be compiled.")
endif()

# Boost and its components
find_package( Boost REQUIRED COMPONENTS program_options )

if ( NOT Boost_FOUND )

  message(STATUS "This project requires the Boost library, and will not be compiled.")

  return()  

endif()

# include for local directory
include_directories( BEFORE include )

# include for local package
include_directories( BEFORE ../../include )

# Boost linking
add_definitions( "-DCGAL_USE_BOOST_PROGRAM_OPTIONS" )
list(APPEND CGAL_3RD_PARTY_LIBRARIES ${Boost_PROGRAM_OPTIONS_LIBRARY} )


# Creating entries for all .cpp/.C files with "main" routine
# ##########################################################

include( CGAL_CreateSingleSourceCGALProgram )

create_single_source_cgal_program( "compute_isolines.cpp" )

create_single_source_cgal_program( "shortest_paths_elephant.cpp" )

create_single_source_cgal_program( "shortest_paths_elephant_with_id.cpp" )

create_single_source_cgal_program( "shortest_paths_elephant_multiple_sources.cpp" )

if(OpenMesh_FOUND)
  create_single_source_cgal_program( "shortest_paths_elephant_OpenMesh.cpp" )
  target_link_libraries( shortest_paths_elephant_OpenMesh ${OPENMESH_LIBRARIES} )
endif()
