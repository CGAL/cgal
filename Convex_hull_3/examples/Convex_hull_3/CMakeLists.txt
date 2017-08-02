# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

project( Convex_hull_3_Examples )


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

find_package( OpenMesh QUIET )

if ( OpenMesh_FOUND )
  include( UseOpenMesh )
else()
  message(STATUS "Examples that use OpenMesh will not be compiled.")
endif()

# include for local directory

# include for local package
include_directories( BEFORE ../../include )


# Creating entries for all C++ files with "main" routine
# ##########################################################

include( CGAL_CreateSingleSourceCGALProgram )

create_single_source_cgal_program( "dynamic_hull_3.cpp" )

create_single_source_cgal_program( "dynamic_hull_SM_3.cpp" )

create_single_source_cgal_program( "halfspace_intersection_3.cpp" )

create_single_source_cgal_program( "lloyd_algorithm.cpp" )

create_single_source_cgal_program( "quickhull_3.cpp" )

create_single_source_cgal_program( "quickhull_any_dim_3.cpp" )


if(OpenMesh_FOUND)
  create_single_source_cgal_program( "quickhull_OM_3.cpp" )
  create_single_source_cgal_program( "dynamic_hull_OM_3.cpp" )

  target_link_libraries( quickhull_OM_3 ${OPENMESH_LIBRARIES} )
  target_link_libraries( dynamic_hull_OM_3 ${OPENMESH_LIBRARIES} )
endif()
