# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.1...3.15)
project( Polygon_mesh_processing )

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

# include for local directory

# include for local package
find_package(Eigen3 3.2.0) #(requires 3.2.0 or greater)
if (EIGEN3_FOUND)
  # Executables that require Eigen 3.2
  include( ${EIGEN3_USE_FILE} )
endif(EIGEN3_FOUND)

# Creating entries for all .cpp/.C files with "main" routine
# ##########################################################

create_single_source_cgal_program( "polygon_mesh_slicer.cpp" )

create_single_source_cgal_program( "hole_filling.cpp" )


