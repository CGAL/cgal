# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.1...3.15)
project( Classification_Tests )


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

set(Classification_dependencies_met TRUE)

find_package( Boost OPTIONAL_COMPONENTS serialization iostreams )
if (NOT Boost_SERIALIZATION_FOUND)
  message(STATUS "NOTICE: This project requires Boost Serialization, and will not be compiled.")
  set(Classification_dependencies_met FALSE)
endif()
if (NOT Boost_IOSTREAMS_FOUND)
  message(STATUS "NOTICE: This project requires Boost IO Streams, and will not be compiled.")
  set(Classification_dependencies_met FALSE)
endif()

find_package(Eigen3 3.1.0 REQUIRED) #(3.1.0 or greater)
if (NOT EIGEN3_FOUND)
  message(STATUS "This project requires the Eigen library, and will not be compiled.")
  set(Classification_dependencies_met FALSE)
endif()

find_package(TBB QUIET)

if (NOT Classification_dependencies_met)
  return()
endif()

create_single_source_cgal_program( "test_classification_point_set.cpp" )
create_single_source_cgal_program( "test_classification_io.cpp" )

foreach(target test_classification_point_set test_classification_io)
  CGAL_target_use_Eigen(${target})
  CGAL_target_use_Boost_IOStreams(${target})
  CGAL_target_use_Boost_Serialization(${target})
  if(TBB_FOUND)
    CGAL_target_use_TBB(${target})
  endif()
endforeach()

