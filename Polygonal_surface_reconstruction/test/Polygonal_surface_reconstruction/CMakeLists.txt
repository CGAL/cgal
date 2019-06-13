# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

project( Polygonal_surface_reconstruction_Tests )


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

find_package( SCIP )


# include for local directory

# include for local package
include_directories( BEFORE ../../include )
include_directories( BEFORE ${SCIP_INCLUDE_DIRS} )


# Creating entries for all C++ files with "main" routine
# ##########################################################

include( CGAL_CreateSingleSourceCGALProgram )

create_single_source_cgal_program( "polygonal_surface_reconstruction_test.cpp" )

if ( SCIP_FOUND )
  target_link_libraries( polygonal_surface_reconstruction_test PRIVATE ${SCIP_LIBRARIES} )
endif()

find_path(GLPK_INCLUDE_DIR
  NAMES glpk.h
  DOC "Path to GLPK headers"
  )
find_library(GLPK_LIBRARIES NAMES GLPK DOC "Path to the GLPK library")

if(GLPK_LIBRARIES AND GLPK_INCLUDE_DIR)
  set( GLPK_FOUND TRUE)
endif()

if ( NOT SCIP_FOUND )
  if ( NOT GLPK_FOUND )
    message(STATUS "This project requires the SCIP library or the GLPK library, and will not be compiled.")
    target_compile_definitions(polygonal_surface_reconstruction_test PRIVATE SKIP_TEST_USING_GLPK)
    message("NOTICE: SCIP and GLPK not found")
  return()  
  endif()
endif()


if ( GLPK_FOUND )
	include_directories( BEFORE ${GLPK_INCLUDE_DIR} )
	target_link_libraries( polygonal_surface_reconstruction_test PRIVATE ${GLPK_LIBRARIES} )
	message("GLPK found and used") 
endif()
