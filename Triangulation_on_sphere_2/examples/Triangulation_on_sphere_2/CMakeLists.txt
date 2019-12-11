# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.15)

project( Triangulation_on_sphere_2_Examples )

find_package(CGAL QUIET)

if ( CGAL_FOUND )

  create_single_source_cgal_program( "constrained_triang_on_sphere_ex.cpp" )
  create_single_source_cgal_program( "triang_on_sphere_ex.cpp" )
  create_single_source_cgal_program( "triang_on_sphere_ex_proj.cpp" )

else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()

