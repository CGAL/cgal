# This is the CMake script for compiling a CGAL application.

project( pdb ) 

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )  
endif()

find_package(CGAL REQUIRED PDB )
include( ${CGAL_USE_FILE} )

INCLUDE_DIRECTORIES ( BEFORE ./include )
set(CMAKE_CXX_FLAGS "-g -Wall")
set(CMAKE_CPP_FLAGS "-I/flute1/home/drussel/fs/include")

if (CGAL_PDB_FOUND )

include(CGAL_CreateSingleSourceCGALProgram)


link_libraries(${CGAL_PDB_LIBRARY})

create_single_source_cgal_program( check_heatom.cpp )
create_single_source_cgal_program( check_pdb.cpp )
create_single_source_cgal_program( check_protein.cpp )
create_single_source_cgal_program( check_transform.cpp )
create_single_source_cgal_program( check_bonds.cpp )
create_single_source_cgal_program( check_refine_alignment.cpp )



else()

  message(STATUS "NOTICE: This demo requires PDB, and will not be compiled.")

endif()
