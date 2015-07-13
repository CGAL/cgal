project( Core_examples )


cmake_minimum_required(VERSION 2.8.11)


set( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true )

if ( COMMAND cmake_policy )

  cmake_policy( SET CMP0003 NEW )

endif()

# CGAL and its components
find_package( CGAL QUIET COMPONENTS Core )

if ( NOT CGAL_Core_FOUND )

  message(STATUS "This project requires the CGAL_Core library, and will not be compiled.")
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

# include for local directory

# include for local package
include_directories( BEFORE ../../include )


# Creating entries for all .cpp/.C files with "main" routine
# ##########################################################

include( CGAL_CreateSingleSourceCGALProgram )

create_single_source_cgal_program( "delaunay.cpp" )
