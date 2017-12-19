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
find_package( Boost REQUIRED COMPONENTS serialization iostreams)

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
  message(STATUS "Found OpenCV ${OpenCV_VERSION}")
  include_directories( ${OpenCV_INCLUDE_DIRS} )
  create_single_source_cgal_program( "example_opencv_random_forest.cpp" CXX_FEATURES ${needed_cxx_features} )
  target_link_libraries( example_opencv_random_forest ${OpenCV_LIBS} )
  target_compile_definitions(example_opencv_random_forest PUBLIC "-DCGAL_LINKED_WITH_OPENCV")
else()
  message(STATUS "OpenCV not found, OpenCV random forest example won't be compiled.")
endif()

if (Boost_SERIALIZATION_FOUND AND Boost_IOSTREAMS_FOUND)
  create_single_source_cgal_program( "example_ethz_random_forest.cpp" CXX_FEATURES ${needed_cxx_features} )
  target_link_libraries( example_ethz_random_forest ${Boost_SERIALIZATION_LIBRARY} ${Boost_IOSTREAMS_LIBRARY})
  target_compile_definitions(example_ethz_random_forest PUBLIC "-DCGAL_LINKED_WITH_BOOST_SERIALIZATION")
else()
  message(STATUS "Boost serialization and IO streams not found, ETHZ random forest example won't be compiled.")
endif()

