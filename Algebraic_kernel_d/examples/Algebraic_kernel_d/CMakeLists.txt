# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

project( Algebraic_kernel_d_examples )


cmake_minimum_required(VERSION 2.8.11)

# CGAL and its components
find_package( CGAL QUIET COMPONENTS Core )

if ( NOT CGAL_FOUND )
  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()  
endif()

# include helper file
include( ${CGAL_USE_FILE} )


# Boost and its components
find_package( Boost REQUIRED )
if ( NOT Boost_FOUND )
  message(STATUS "This project requires the Boost library, and will not be compiled.")
  return()
endif()

find_package( MPFI QUIET )
if( MPFI_FOUND )
  include( ${MPFI_USE_FILE} )
endif()

include_directories( BEFORE ../../include )

include( CGAL_CreateSingleSourceCGALProgram )


create_single_source_cgal_program( "Compare_1.cpp" )
create_single_source_cgal_program( "Construct_algebraic_real_1.cpp" )
create_single_source_cgal_program( "Isolate_1.cpp" )
create_single_source_cgal_program( "Sign_at_1.cpp" )
create_single_source_cgal_program( "Solve_1.cpp" )


find_package( RS )
find_package( RS3 )
if ( NOT RS3_FOUND )
  message(STATUS "Parts of this project require the RS libraries, and will not be compiled.")
else()
  include( ${RS_USE_FILE} )
  include( ${RS3_USE_FILE} )
  create_single_source_cgal_program( "RS_Compare_1.cpp" )
  create_single_source_cgal_program( "RS_Construct_algebraic_real_1.cpp" )
  create_single_source_cgal_program( "RS_Isolate_1.cpp" )
  create_single_source_cgal_program( "RS_Sign_at_1.cpp" )
  create_single_source_cgal_program( "RS_Solve_1.cpp" )
endif()
