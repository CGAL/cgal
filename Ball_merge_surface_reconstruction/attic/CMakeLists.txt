# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.23)
project( attic )

find_package(CGAL REQUIRED QUIET OPTIONAL_COMPONENTS Core )

include_directories (BEFORE "../include")

create_single_source_cgal_program( "AutomaticBallMerge.cpp" )
create_single_source_cgal_program( "BallMergeParallel.cpp" )

