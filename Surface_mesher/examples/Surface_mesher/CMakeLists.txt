# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project( Surface_mesher_example ) 

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

find_package(CGAL REQUIRED ImageIO )


include( ${CGAL_USE_FILE} )
include( CreateSingleSourceCGALProgram )

create_single_source_cgal_program( "mesh_a_3d_gray_image.cpp" )
create_single_source_cgal_program( "mesh_an_implicit_function.cpp" )
create_single_source_cgal_program( "polyhedron_remesher.cpp" )

