# Created by the script cgal_create_cmake_script_with_options
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.1...3.15)
project( Triangulation_apps )

# CGAL and its components
find_package( CGAL QUIET COMPONENTS  )

if ( NOT CGAL_FOUND )

  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()  

endif()

# Boost and its components
find_package( Boost REQUIRED )

if ( NOT Boost_FOUND )

  message(STATUS "This project requires the Boost library, and will not be compiled.")

  return()  

endif()

find_package(Eigen3 3.1.0) 
if (EIGEN3_FOUND) 
  include( ${EIGEN3_USE_FILE} ) 
endif()

# include for local directory
include_directories( BEFORE include )

# include for local package

# Creating entries for all .cpp/.C files with "main" routine
# ##########################################################

create_single_source_cgal_program( "points_to_RT_to_off.cpp" )
create_single_source_cgal_program( "points_to_DT_to_off.cpp" )


