# This is the CMake script for compiling a CGAL application.

project( pdb ) 

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

find_package(CGAL REQUIRED PDB )
include( ${CGAL_USE_FILE} )


if (CGAL_FOUND)

include(CreateSingleSourceCGALProgram)

create_single_source_cgal_program( pdb_cat.cpp )
create_single_source_cgal_program( pdb_echo.cpp )
create_single_source_cgal_program( pdb_split.cpp )
create_single_source_cgal_program( pdb_transform.cpp )
create_single_source_cgal_program( pdb_tree.cpp )


else()

  message(STATUS "NOTICE:his demo requires CGAL and PDB, and will not be compiled.")

endif()
