# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

project( Surface_mesh_approximation )


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

# include for local directory

# include for local package
include_directories( BEFORE ../../include )


# Creating entries for all C++ files with "main" routine
# ##########################################################

include( CGAL_CreateSingleSourceCGALProgram )

create_single_source_cgal_program( "vsa_class_interface_test.cpp" )

create_single_source_cgal_program( "vsa_correctness_test.cpp" )

create_single_source_cgal_program( "vsa_error_decrease_test.cpp" )

create_single_source_cgal_program( "vsa_kernel_test.cpp" )

create_single_source_cgal_program( "vsa_approximation_test.cpp" )

create_single_source_cgal_program( "vsa_segmentation_test.cpp" )

create_single_source_cgal_program( "vsa_meshing_manifold_test.cpp" )

create_single_source_cgal_program( "vsa_metric_test.cpp" )

create_single_source_cgal_program( "vsa_teleportation_test.cpp" )


