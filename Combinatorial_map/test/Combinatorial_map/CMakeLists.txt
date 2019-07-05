# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Combinatorial_map_Tests )

cmake_minimum_required(VERSION 3.1)

find_package(CGAL QUIET)

if ( CGAL_FOUND )

  # create a target per cppfile
  file(GLOB cppfiles RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
  foreach(cppfile ${cppfiles})
    create_single_source_cgal_program( "${cppfile}" )
  endforeach()

  find_package( OpenMesh QUIET )
  if (TARGET OpenMesh::OpenMesh)
    target_link_libraries(Combinatorial_map_copy_test PRIVATE OpenMesh::OpenMesh)
  endif()

else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()

