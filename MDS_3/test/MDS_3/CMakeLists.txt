# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.15)
project( MDS_3_Tests )

if(POLICY CMP0074)
  cmake_policy(SET CMP0074 NEW)
endif()

find_package(CGAL REQUIRED)

# Use Eigen
find_package(Eigen3 3.1.0 REQUIRED) #(3.1.0 or greater)
include(CGAL_Eigen3_support)
if (NOT TARGET CGAL::Eigen3_support)
  message(STATUS "This project requires the Eigen library, and will not be compiled.")
  return()
endif()

create_single_source_cgal_program( "test_c3t3.cpp" )
create_single_source_cgal_program( "test_c3t3_io.cpp" )
create_single_source_cgal_program( "test_c3t3_with_features.cpp" )
create_single_source_cgal_program( "test_c3t3_into_facegraph.cpp" )
create_single_source_cgal_program( "test_c3t3_extract_subdomains_boundaries.cpp" )

foreach(target
    test_c3t3
    test_c3t3_io
    test_c3t3_with_features
    test_c3t3_into_facegraph
    test_c3t3_extract_subdomains_boundaries)
  if(TARGET ${target})
    target_link_libraries(${target} PUBLIC CGAL::Eigen3_support)
  endif()
endforeach()
