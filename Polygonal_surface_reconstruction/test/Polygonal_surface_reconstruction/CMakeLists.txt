# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

project( Polygonal_surface_reconstruction_Tests )


cmake_minimum_required(VERSION 3.1...3.15)

# CGAL and its components
find_package( CGAL QUIET COMPONENTS  )

if ( NOT CGAL_FOUND )

  message(STATUS "NOTICE: This project requires the CGAL library, and will not be compiled.")
  return()  

endif()

# include helper file
include( ${CGAL_USE_FILE} )


# Boost and its components
find_package( Boost REQUIRED )

if ( NOT Boost_FOUND )

  message(STATUS "NOTICE: This project requires the Boost library, and will not be compiled.")

  return()  

endif()

# Creating entries for all C++ files with "main" routine
# ##########################################################

include( CGAL_CreateSingleSourceCGALProgram )


find_package(Eigen3 3.1.0) #(requires 3.1.0 or greater)
if(NOT EIGEN3_FOUND)
  message(STATUS "NOTICE: This project requires Eigen 3.1 (or greater) and will not be compiled.")
  return()
endif()

find_package(SCIP QUIET)
if (NOT SCIP_FOUND )
  find_package( GLPK QUIET)
  if ( NOT GLPK_FOUND )
    message( STATUS "NOTICE : This project requires either SCIP or GLPK, and will not be compiled.")
    return()
  endif()
endif()

create_single_source_cgal_program( "polygonal_surface_reconstruction_test.cpp")
CGAL_target_use_Eigen(polygonal_surface_reconstruction_test)
if (SCIP_FOUND)
  CGAL_target_use_SCIP(polygonal_surface_reconstruction_test)
else()
  CGAL_target_use_GLPK(polygonal_surface_reconstruction_test)
endif()
