# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project( Periodic_3_mesh_3_Tests )

cmake_minimum_required( VERSION 2.8.11 )

find_package( CGAL QUIET )

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )

  include_directories( BEFORE ../../include )

  # Use Eigen
  find_package(Eigen3 3.1.0 REQUIRED) #(3.1.0 or greater)
  if (NOT EIGEN3_FOUND)
    message(STATUS "This project requires the Eigen library, and will not be compiled.")
    return()
  else()
    include( ${EIGEN3_USE_FILE} )
  endif()

  include( CGAL_CreateSingleSourceCGALProgram )
  create_single_source_cgal_program( "test_implicit_shapes_bunch.cpp" )
  create_single_source_cgal_program( "test_implicit_shapes_with_features.cpp" )
  create_single_source_cgal_program( "test_triply_periodic_minimal_surfaces.cpp" )

else()

  message(STATUS "This program requires the CGAL library, and will not be compiled.")

endif()

