# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Surface_mesh_ )

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

  include( CGAL_CreateSingleSourceCGALProgram )

  include_directories (BEFORE "../../include")

  create_single_source_cgal_program( "sm_join.cpp" )
  create_single_source_cgal_program( "sm_aabbtree.cpp" )
  create_single_source_cgal_program( "sm_bgl.cpp" )
  create_single_source_cgal_program( "sm_circulators.cpp" )
  create_single_source_cgal_program( "sm_do_intersect.cpp" )
  create_single_source_cgal_program( "sm_iterators.cpp" )
  create_single_source_cgal_program( "sm_kruskal.cpp" )
  create_single_source_cgal_program( "sm_memory.cpp" )
  create_single_source_cgal_program( "sm_properties.cpp" )

else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()

