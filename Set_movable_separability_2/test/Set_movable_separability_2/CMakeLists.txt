# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Set_movable_separability_2_Tests )

cmake_minimum_required(VERSION 3.1)

# Use C++11 for this directory and its sub-directories.
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)

find_package(CGAL QUIET COMPONENTS Core )

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )

  include( CGAL_CreateSingleSourceCGALProgram )

  include_directories (BEFORE "../../include")

create_single_source_cgal_program( "test_top_edges_single_mold_translational_casting.cpp" )
create_single_source_cgal_program( "test_is_pullout_directions_single_mold_translational_casting.cpp" )

else()

    message(STATUS "This program requires the CGAL library, and will not be compiled.")

endif()
