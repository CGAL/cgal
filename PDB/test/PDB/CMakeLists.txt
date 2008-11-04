# This is the CMake script for compiling a CGAL application.

project( pdb ) 

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )  
endif()

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

find_package(CGAL REQUIRED PDB )
include( ${CGAL_USE_FILE} )

INCLUDE_DIRECTORIES ( BEFORE ./include )

if (CGAL_PDB_FOUND )

include(CreateSingleSourceCGALProgram)

create_single_source_cgal_program( check_bonds.cpp )
create_single_source_cgal_program( check_heatom.cpp )
create_single_source_cgal_program( check_pdb.cpp )
create_single_source_cgal_program( check_protein.cpp )
create_single_source_cgal_program( check_refine_alignment.cpp )
create_single_source_cgal_program( check_transform.cpp )
create_single_source_cgal_program( check_transforms.cpp )



else()

  message(STATUS "NOTICE: This demo requires PDB, and will not be compiled.")

endif()
