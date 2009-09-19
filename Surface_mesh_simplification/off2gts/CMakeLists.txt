# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Off2gts )

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)
 
if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )  
endif()
 
find_package(CGAL REQUIRED )

include( ${CGAL_USE_FILE} )

include( CGAL_CreateSingleSourceCGALProgram )

create_single_source_cgal_program( "off2gts.cpp" )


