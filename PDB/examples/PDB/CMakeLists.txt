# This is the CMake script for compiling a CGAL application.

project( PDB_example ) 

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

find_package(CGAL REQUIRED PDB )
include( ${CGAL_USE_FILE} )

include( CreateSingleSourceCGALProgram )

create_single_source_cgal_program( "extracting_geometry.cpp" )

