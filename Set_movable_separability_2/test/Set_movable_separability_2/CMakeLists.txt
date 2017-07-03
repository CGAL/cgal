# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Casting_2_test )

cmake_minimum_required(VERSION 2.8.10)

find_package(CGAL QUIET COMPONENTS Core )

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )

  include( CGAL_CreateSingleSourceCGALProgram )

  include_directories (BEFORE "../../include")

  if (CMAKE_COMPILER_IS_GNUCXX)
    add_definitions(-std=c++11)
  endif()

create_single_source_cgal_program( "test_top_edges_single_mold_translational_casting.cpp" )
create_single_source_cgal_program( "test_is_pullout_directions_single_mold_translational_casting.cpp" )

else()

    message(STATUS "This program requires the CGAL library, and will not be compiled.")

endif()
