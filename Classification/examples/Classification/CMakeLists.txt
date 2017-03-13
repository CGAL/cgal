# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

project( Classification )


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

find_package( TBB REQUIRED )
if( TBB_FOUND )
  include(${TBB_USE_FILE})
  list(APPEND CGAL_3RD_PARTY_LIBRARIES ${TBB_LIBRARIES})
endif()


# include for local directory
include_directories( BEFORE include )

# include for local package
include_directories( BEFORE ../../include )


# Creating entries for all C++ files with "main" routine
# ##########################################################

include( CGAL_CreateSingleSourceCGALProgram )

create_single_source_cgal_program( "example_classifier.cpp" )
create_single_source_cgal_program( "example_point_set_classifier.cpp" )
create_single_source_cgal_program( "example_feature.cpp" )
