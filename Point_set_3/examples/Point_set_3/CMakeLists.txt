# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

project( Point_set_3_Examples )


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
include_directories( BEFORE include )

# include for local package
include_directories( BEFORE ../../include )


# Creating entries for all C++ files with "main" routine
# ##########################################################

include( CGAL_CreateSingleSourceCGALProgram )

create_single_source_cgal_program( "point_set.cpp" )
create_single_source_cgal_program( "point_set_property.cpp" )
create_single_source_cgal_program( "point_set_read_xyz.cpp" )
create_single_source_cgal_program( "point_set_advanced.cpp" )

set(needed_cxx_features cxx_rvalue_references cxx_variadic_templates)    
create_single_source_cgal_program( "point_set_read_ply.cpp" CXX_FEATURES ${needed_cxx_features} )

find_package(Eigen3 3.1.0) #(requires 3.1.0 or greater)
if (NOT EIGEN3_FOUND)
  find_package(LAPACK)
  if(LAPACK_FOUND)
    include( ${LAPACK_USE_FILE} )
  endif(LAPACK_FOUND)
else()
  include( ${EIGEN3_USE_FILE} )
  create_single_source_cgal_program( "point_set_algo.cpp" )
endif()



