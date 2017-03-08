# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project( Periodic_3_mesh_3_Tests )

cmake_minimum_required(VERSION 2.8.10)

set(CMAKE_CXX_STANDARD 14)

find_package(CGAL QUIET COMPONENTS Core )

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )

  include( CGAL_CreateSingleSourceCGALProgram )

  include_directories (BEFORE "../../include")

#  create_single_source_cgal_program( "combined_spheres.cpp" )
#  create_single_source_cgal_program( "implicit_surface_mesher_test.cpp" )
#  create_single_source_cgal_program( "mesh_moebius_strip.cpp" )
#  create_single_source_cgal_program( "test_c2t3_iterators.cpp" )
#  create_single_source_cgal_program( "test_robust_circumcenter.cpp" )

else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()

