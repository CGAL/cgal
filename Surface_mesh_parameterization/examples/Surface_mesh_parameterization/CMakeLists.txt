# This is the CMake script for compiling this folder.


project( Surface_mesh_parameterization_Examples )

cmake_minimum_required(VERSION 2.8.11)


# Include this package's headers first
include_directories (BEFORE . include ../../include)

# Find CGAL
find_package(CGAL QUIET)

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )
  include( CGAL_CreateSingleSourceCGALProgram )
  find_package(Eigen3 3.1.0) #(requires 3.1.0 or greater)
 
  if (EIGEN3_FOUND)
    # Executables that require Eigen 3.1
    include( ${EIGEN3_USE_FILE} ) 

  create_single_source_cgal_program( "lscm.cpp" )
  create_single_source_cgal_program( "seam.cpp" )
  create_single_source_cgal_program( "discrete_authalic.cpp" )
  create_single_source_cgal_program( "Simple_parameterization.cpp" )

  else(EIGEN3_FOUND)
    message(STATUS "NOTICE: The examples require Eigen 3.1 (or greater) and will not be compiled.")  
  endif(EIGEN3_FOUND)

else()

  message(STATUS "NOTICE: This program requires the CGAL library, and will not be compiled.")

endif()

