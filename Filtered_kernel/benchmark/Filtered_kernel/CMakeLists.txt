# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Filtered_kernel_test )

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)
 
if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )  
endif()
 
find_package(CGAL QUIET COMPONENTS Core )

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )

  include( CGAL_CreateSingleSourceCGALProgram )

  include_directories (BEFORE ../../include)

  option(ONLY_TEST_COMPARISONS "Check comparisons of doubles instead of Orientation_3")

  if(ONLY_TEST_COMPARISONS)
    add_definitions(-DONLY_TEST_COMPARISONS)
  endif()
  create_single_source_cgal_program( "orientation_3.cpp" )


else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()

