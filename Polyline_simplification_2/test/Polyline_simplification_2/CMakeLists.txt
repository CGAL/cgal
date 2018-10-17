# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

project( Polyline_simplification_2_Tests )


cmake_minimum_required(VERSION 2.8.11)



# CGAL and its components
find_package( CGAL QUIET COMPONENTS  )

if ( NOT CGAL_FOUND )

  message(STATUS "This project requires the CGAL library, and will not be compiled.")
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

# Creating entries for all .cpp/.C files with "main" routine
# ##########################################################

create_single_source_cgal_program( "simplify_polygon_test.cpp" )

create_single_source_cgal_program( "simplify_polyline_with_duplicate_points.cpp" )


