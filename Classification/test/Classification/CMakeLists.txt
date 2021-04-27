# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.1...3.20)
project(Classification_Tests)

# CGAL and its components
find_package(CGAL REQUIRED)

# Boost and its components
find_package(Boost REQUIRED)

if(NOT Boost_FOUND)
  message(
    STATUS "This project requires the Boost library, and will not be compiled.")
  return()
endif()

set(Classification_dependencies_met TRUE)

find_package(Boost OPTIONAL_COMPONENTS serialization iostreams)
include(CGAL_Boost_serialization_support)
include(CGAL_Boost_iostreams_support)

if(NOT TARGET CGAL::Boost_serialization_support)
  message(
    STATUS
      "NOTICE: This project requires Boost Serialization, and will not be compiled."
  )
  set(Classification_dependencies_met FALSE)
endif()
if(NOT TARGET CGAL::Boost_iostreams_support)
  message(
    STATUS
      "NOTICE: This project requires Boost IO Streams, and will not be compiled."
  )
  set(Classification_dependencies_met FALSE)
endif()

find_package(Eigen3 3.1.0 REQUIRED) #(3.1.0 or greater)
include(CGAL_Eigen3_support)
if(NOT TARGET CGAL::Eigen3_support)
  message(
    STATUS "This project requires the Eigen library, and will not be compiled.")
  set(Classification_dependencies_met FALSE)
endif()

find_package(TBB QUIET)
include(CGAL_TBB_support)

if(NOT Classification_dependencies_met)
  return()
endif()

create_single_source_cgal_program("test_classification_point_set.cpp")
create_single_source_cgal_program("test_classification_io.cpp")

foreach(target test_classification_point_set test_classification_io)
  target_link_libraries(
    ${target} PUBLIC CGAL::Eigen3_support CGAL::Boost_iostreams_support
                     CGAL::Boost_serialization_support)
  if(TARGET CGAL::TBB_support)
    target_link_libraries(${target} PUBLIC CGAL::TBB_support)
  endif()
endforeach()
