project( Core_Examples )


cmake_minimum_required(VERSION 2.8.11)



# CGAL and its components
find_package( CGAL QUIET COMPONENTS Core )

# include helper file
include( ${CGAL_USE_FILE} )

if ( NOT CGAL_Core_FOUND )

  message(STATUS "This project requires the CGAL_Core library, and will not be compiled.")
  return()

endif()

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
