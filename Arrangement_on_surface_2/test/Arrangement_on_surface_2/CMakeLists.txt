# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Arrangement_on_surface_2_Tests )

enable_testing()

cmake_minimum_required(VERSION 3.1)

find_package(CGAL QUIET COMPONENTS Core)

if ( CGAL_FOUND )

  if(COMMAND target_compile_options)
    # Since CMake-2.8.12: New CMake script, that defines the targets and
    # the CTest test cases.
    include( ${CMAKE_CURRENT_SOURCE_DIR}/cgal_test.cmake )
  else()
    # If CMake version is <= 2.8.11, use the usual CMake script.

    # create a target per cppfile
    file(GLOB cppfiles RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
    foreach(cppfile ${cppfiles})
      create_single_source_cgal_program( "${cppfile}" )
    endforeach()
  endif()

else()

    message(STATUS "This program requires the CGAL library, and will not be compiled.")

endif()
