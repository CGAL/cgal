# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.1...3.14)

project( Tetrahedral_remeshing_Examples )

# CGAL and its components
find_package( CGAL REQUIRED )
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


# Creating entries for all C++ files with "main" routine
# ##########################################################
create_single_source_cgal_program( "tetrahedral_remeshing_example.cpp" )
create_single_source_cgal_program( "tetrahedral_remeshing_with_features.cpp")
create_single_source_cgal_program( "tetrahedral_remeshing_of_one_subdomain.cpp")
