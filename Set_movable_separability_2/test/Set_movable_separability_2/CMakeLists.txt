# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project( Set_movable_separability_2_Tests )

cmake_minimum_required(VERSION 3.1)

if (RUNNING_CGAL_AUTO_TEST)
  # Just to avoid a warning from CMake when that variable is set on the command line...
endif()
if (CGAL_DIR)
  # Just to avoid a warning from CMake when that variable is set on the command line...
endif()

list(FIND CMAKE_CXX_COMPILE_FEATURES cxx_generalized_initializers has_cpp11)
if (has_cpp11 LESS 0)
  message(STATUS "NOTICE: These examples requires a C++11 compiler and will not be compiled.")
  return()
endif()

# Use C++11 for this directory and its sub-directories.
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)

find_package(CGAL QUIET COMPONENTS )

if (CGAL_FOUND)
  include( ${CGAL_USE_FILE} )
  include( CGAL_CreateSingleSourceCGALProgram )
  include_directories (BEFORE "../../include")
  create_single_source_cgal_program( "test_top_edges_single_mold_trans_cast.cpp" )
  create_single_source_cgal_program( "test_is_pullout_directions_single_mold_trans_cast.cpp" )

else()
  message(STATUS "This program requires the CGAL library, and will not be compiled.")
endif()
