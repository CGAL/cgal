# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.
# Then modified by hand to add Eigen3.


cmake_minimum_required(VERSION 3.1...3.15)
project( NewKernel_d_Tests )


if(CMAKE_COMPILER_IS_GNUCCX AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.4)
  message(STATUS "NOTICE: this directory requires a version of gcc >= 4.4, and will not be compiled.")
  return()
endif()

find_package(CGAL QUIET)

if ( CGAL_FOUND )

  find_package(Eigen3)
  if (EIGEN3_FOUND)
    create_single_source_cgal_program( "Epick_d.cpp" )
    CGAL_target_use_Eigen(Epick_d)
  else()

    message(STATUS "NOTICE: This program requires the Eigen3 library, and will not be compiled.")

  endif()

else()

    message(STATUS "NOTICE: This program requires the CGAL library, and will not be compiled.")

endif()

