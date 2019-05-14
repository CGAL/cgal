# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

project( Periodic_3_mesh_3_Examples )

cmake_minimum_required( VERSION 2.8.11 )

# CGAL and its components
find_package( CGAL QUIET )

if ( NOT CGAL_FOUND )
  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()
endif()

# include helper file
include( ${CGAL_USE_FILE} )

# include for local package
include_directories( BEFORE ../../include )

# Use Eigen
find_package(Eigen3 3.1.0 REQUIRED) #(3.1.0 or greater)
if (NOT EIGEN3_FOUND)
  message(STATUS "This project requires the Eigen library, and will not be compiled.")
  return()
else()
  include( ${EIGEN3_USE_FILE} )
endif()

# Boost and its components
find_package( Boost )

if ( NOT Boost_FOUND )
  message(STATUS "This project requires the Boost library, and will not be compiled.")
  return()
endif()

# Creating entries for all C++ files with "main" routine
# ##########################################################

include( CGAL_CreateSingleSourceCGALProgram )

create_single_source_cgal_program( "mesh_implicit_shape.cpp" )
create_single_source_cgal_program( "mesh_implicit_multi_domain.cpp" )
create_single_source_cgal_program( "mesh_implicit_shape_with_subdomains.cpp" )
create_single_source_cgal_program( "mesh_implicit_shape_with_optimizers.cpp" )
create_single_source_cgal_program( "mesh_implicit_shape_with_features.cpp" )
