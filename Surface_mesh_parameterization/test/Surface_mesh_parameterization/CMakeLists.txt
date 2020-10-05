# This is the CMake script for compiling this folder.


cmake_minimum_required(VERSION 3.1...3.15)
project( Surface_mesh_parameterization_Tests )



# Find CGAL
find_package(CGAL QUIET)

if ( CGAL_FOUND )


  find_package(Eigen3 3.1.0) #(requires 3.1.0 or greater)
  include(CGAL_Eigen_support)
  if(TARGET CGAL::Eigen_support)
    create_single_source_cgal_program( "extensive_parameterization_test.cpp" )
    target_link_libraries(extensive_parameterization_test PUBLIC CGAL::Eigen_support)
  else()
    message(STATUS "NOTICE: The tests require Eigen 3.1 (or greater) and will not be compiled.")
  endif()

else()

  message(STATUS "NOTICE: This program requires the CGAL library, and will not be compiled.")

endif()
