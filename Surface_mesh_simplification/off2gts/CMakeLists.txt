# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Off2gts )

CMAKE_MINIMUM_REQUIRED(VERSION 2.6.2)

find_package(CGAL REQUIRED )

include( ${CGAL_USE_FILE} )

include( CGAL_CreateSingleSourceCGALProgram )

create_single_source_cgal_program( "off2gts.cpp" )


