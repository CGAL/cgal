# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project( Skin_surface_3_example ) 

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )  
endif()

find_package(CGAL QUIET COMPONENTS PDB )
include( ${CGAL_USE_FILE} )

if(CGAL_PDB_FOUND)

  include( CGAL_CreateSingleSourceCGALProgram )

  include_directories (BEFORE include)

  create_single_source_cgal_program( "NGHK_skin_surface_simple.cpp" )
  create_single_source_cgal_program( "NGHK_skin_surface_subdiv.cpp" )
  create_single_source_cgal_program( "skin_surface_pdb_reader.cpp" )
  create_single_source_cgal_program( "skin_surface_simple.cpp" )
  create_single_source_cgal_program( "skin_surface_subdiv.cpp" )
  create_single_source_cgal_program( "skin_surface_subdiv_with_normals.cpp" )
  create_single_source_cgal_program( "union_of_balls_simple.cpp" )
  create_single_source_cgal_program( "union_of_balls_subdiv.cpp" )

else(CGAL_PDB_FOUND)

  message(STATUS "NOTICE: This directory requires CGAL PDB library, and will not be compiled.")

endif(CGAL_PDB_FOUND)
