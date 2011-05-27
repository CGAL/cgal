# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

# cmake ../ -DCMAKE_BUILD_TYPE=Debug 
# ou
# cmake ../ -DCMAKE_BUILD_TYPE=Release

project( Map_3_examples )

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -w -Wall -pedantic") # -Wextra")

if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )
endif()

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0")
endif()

# Pour le problème de valgrind avec CGAL
# add_definition(-DCGAL_DISABLE_ROUNDING_MATH_CHECK)

find_package(CGAL QUIET COMPONENTS Core )

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )

  include( CGAL_CreateSingleSourceCGALProgram )
  include_directories(BEFORE ../../include)

  create_single_source_cgal_program( "map_3_simple_example.cpp" )
  create_single_source_cgal_program( "map_4_simple_example.cpp" )
  create_single_source_cgal_program( "map_3_with_colored_facets.cpp" )
  create_single_source_cgal_program( "map_3_operations.cpp" )
  create_single_source_cgal_program( "map_3_marks.cpp" )
  create_single_source_cgal_program( "test-tuple.cc" )
else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()

