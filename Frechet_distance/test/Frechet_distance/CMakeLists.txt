# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.14)
project( Frechet_distance_Tests )

find_package(CGAL REQUIRED)

find_package(Eigen3 3.1.0 QUIET) #(3.1.0 or greater)
include(CGAL_Eigen3_support)
if(TARGET CGAL::Eigen3_support)
  create_single_source_cgal_program( "Frechet_distance_test.cpp" )
  target_link_libraries(Frechet_distance_test PRIVATE CGAL::Eigen3_support)
  create_single_source_cgal_program( "Frechet_distance_minimal_traits_2.cpp" )
  target_link_libraries(Frechet_distance_minimal_traits_2 PRIVATE CGAL::Eigen3_support)
else()
  message("NOTICE: All tests require the Eigen3 library, and will not be compiled.")
endif()


