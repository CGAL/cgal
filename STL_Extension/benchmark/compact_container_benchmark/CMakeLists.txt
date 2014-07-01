# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project( Compact_container_benchmark )

cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" VERSION_GREATER 2.6)
  if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
    cmake_policy(VERSION 2.8.4)
  else()
    cmake_policy(VERSION 2.6)
  endif()
endif()

find_package(CGAL QUIET COMPONENTS Core )

if ( CGAL_FOUND )
  include( ${CGAL_USE_FILE} )
    
  find_package( TBB )
  
  if( TBB_FOUND )
    include(${TBB_USE_FILE})
    list(APPEND CGAL_3RD_PARTY_LIBRARIES ${TBB_LIBRARIES})
  endif()
  
  include( CGAL_CreateSingleSourceCGALProgram )

  include_directories (BEFORE "../../include")

  create_single_source_cgal_program( "cc_benchmark.cpp" )

else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()
