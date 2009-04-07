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

if (CGAL_PDB_FOUND )

include(CGAL_CreateSingleSourceCGALProgram)


link_libraries(${CGAL_PDB_LIBRARY})

create_single_source_cgal_program(check_atom_hetatom.cpp)
create_single_source_cgal_program(check_bonds.cpp)
create_single_source_cgal_program(check_hetatom.cpp)
create_single_source_cgal_program(check_iterators.cpp)
create_single_source_cgal_program(check_pdb.cpp)
create_single_source_cgal_program(check_refine_alignment.cpp)
create_single_source_cgal_program(check_small_map.cpp)
create_single_source_cgal_program(check_transform.cpp)



else()

  message(STATUS "NOTICE: This demo requires PDB, and will not be compiled.")

endif()
