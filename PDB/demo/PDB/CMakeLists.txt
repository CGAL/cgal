# This is the CMake script for compiling a CGAL application.

project( pdb ) 

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )  
endif()

find_package(CGAL COMPONENTS PDB )
include( ${CGAL_USE_FILE} )

find_package ( Boost COMPONENTS program_options )

if (CGAL_PDB_FOUND  AND  Boost_PROGRAM_OPTIONS_FOUND )

  link_libraries(${Boost_LIBRARIES})

  include(CGAL_CreateSingleSourceCGALProgram)

  create_single_source_cgal_program( pdb_cat.cpp )
  create_single_source_cgal_program( pdb_echo.cpp )
  create_single_source_cgal_program( pdb_split.cpp )
  create_single_source_cgal_program( pdb_transform.cpp )
  create_single_source_cgal_program( pdb_tree.cpp )
  create_single_source_cgal_program( pdb_info.cpp )

else()

  message(STATUS "NOTICE: This demo requires PDB and Boost program_options, and will not be compiled.")

endif()
