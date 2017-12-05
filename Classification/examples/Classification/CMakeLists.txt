# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

project( Classification_Examples )


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

find_package( TBB )

find_package(OpenCV QUIET)

# include for local directory
include_directories( BEFORE include )

# include for local package
include_directories( BEFORE ../../include )


# Creating entries for all C++ files with "main" routine
# ##########################################################

include( CGAL_CreateSingleSourceCGALProgram )

set(needed_cxx_features cxx_rvalue_references cxx_variadic_templates)
create_single_source_cgal_program( "example_classification.cpp" CXX_FEATURES ${needed_cxx_features} )
if(TBB_FOUND)
  CGAL_target_use_TBB(example_classification)
endif()
create_single_source_cgal_program( "example_generation_and_training.cpp" CXX_FEATURES ${needed_cxx_features} )
create_single_source_cgal_program( "example_feature.cpp" CXX_FEATURES ${needed_cxx_features} )

if( OpenCV_FOUND )
  include_directories( ${OpenCV_INCLUDE_DIRS} )
  create_single_source_cgal_program( "example_random_forest.cpp" CXX_FEATURES ${needed_cxx_features} )
  target_link_libraries( example_random_forest PRIVATE ${OpenCV_LIBS} )
else()
  message(STATUS "OpenCV not found, random forest example won't be compiled.")
endif()

