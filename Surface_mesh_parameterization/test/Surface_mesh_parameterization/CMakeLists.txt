# This is the CMake script for compiling this folder.


project( Surface_mesh_parameterization_Tests )

cmake_minimum_required(VERSION 3.1)


# Find CGAL
find_package(CGAL QUIET)

if ( CGAL_FOUND )


  find_package(Eigen3 3.1.0) #(requires 3.1.0 or greater)
  if(EIGEN3_FOUND)
    include( ${EIGEN3_USE_FILE} )
    create_single_source_cgal_program( "extensive_parameterization_test.cpp" )
  else(EIGEN3_FOUND)
    message(STATUS "NOTICE: The tests require Eigen 3.1 (or greater) and will not be compiled.")
  endif(EIGEN3_FOUND)

else()

  message(STATUS "NOTICE: This program requires the CGAL library, and will not be compiled.")

endif()

