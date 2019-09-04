# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.1...3.15)
project( Heat_method_3_Examples )



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


find_package(Eigen3 3.3.0)

if (EIGEN3_FOUND)
  include( ${EIGEN3_USE_FILE} )
else()
  message(STATUS "This project requires the Eigen library (3.3 or greater), and will not be compiled.")

  return()
endif()


# include for local directory
include_directories( BEFORE include )



# Creating entries for all C++ files with "main" routine
# ##########################################################

include( CGAL_CreateSingleSourceCGALProgram )

create_single_source_cgal_program( "heat_method.cpp" )
create_single_source_cgal_program( "heat_method_polyhedron.cpp" )
create_single_source_cgal_program( "heat_method_surface_mesh.cpp" )
create_single_source_cgal_program( "heat_method_surface_mesh_direct.cpp" )
