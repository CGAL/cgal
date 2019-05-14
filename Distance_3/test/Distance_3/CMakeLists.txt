# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Distance_3_Tests )

cmake_minimum_required(VERSION 3.1)

find_package(CGAL QUIET)

if ( CGAL_FOUND )

 if (MSVC)
    # Turn off a VC++ warning on a potential division by zero
    # in Cartesian_kernel/include/CGAL/constructions/kernel_ftC3.h
    # where CGAL_assume() does not help
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4723")
  endif()
  
  # create a target per cppfile
  file(GLOB cppfiles RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
  foreach(cppfile ${cppfiles})
    create_single_source_cgal_program( "${cppfile}" )
  endforeach()

else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()

