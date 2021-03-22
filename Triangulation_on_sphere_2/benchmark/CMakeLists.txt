# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project( Triangulation_on_sphere_2_Benchmarks )

cmake_minimum_required(VERSION 3.1...3.15)

find_package(CGAL REQUIRED COMPONENTS Core )

if ( CGAL_FOUND )

  create_single_source_cgal_program( "bench_dtos2.cpp" )
  create_single_source_cgal_program( "generate_points.cpp" )

else()

    message(STATUS "This program requires the CGAL library, and will not be compiled.")

endif()

