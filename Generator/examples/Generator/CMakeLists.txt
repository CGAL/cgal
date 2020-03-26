# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


cmake_minimum_required(VERSION 3.1...3.15)
project( Generator_Examples )


find_package(CGAL QUIET)

if ( CGAL_FOUND )

  # Use Eigen
  find_package(Eigen3 3.1.0 QUIET) #(3.1.0 or greater)
  if (EIGEN3_FOUND)
    include( ${EIGEN3_USE_FILE} )
  endif()

  # create a target per cppfile
  file(GLOB cppfiles RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
  foreach(cppfile ${cppfiles})
    if(NOT (${cppfile} STREQUAL "random_points_in_tetrahedral_mesh_3.cpp")
        OR NOT (${cppfile} STREQUAL "random_points_on_tetrahedral_mesh_3.cpp")
        OR EIGEN3_FOUND)
      create_single_source_cgal_program( "${cppfile}" )
    endif()
  endforeach()

else()

    message(STATUS "This program requires the CGAL library, and will not be compiled.")

endif()

