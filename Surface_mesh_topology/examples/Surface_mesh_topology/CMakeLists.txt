# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project(Surface_mesh_topology_Examples)

cmake_minimum_required(VERSION 3.1)

find_package(CGAL QUIET COMPONENTS Qt5)
include("CMakeBasicViewerQt.inc")

if(CGAL_Qt5_FOUND)
  add_definitions(-DCGAL_USE_BASIC_VIEWER -DQT_NO_KEYWORDS)
endif()

if ( CGAL_FOUND )

  # create a target per cppfile
 # file(GLOB cppfiles RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
 # foreach(cppfile ${cppfiles})
 #   create_single_source_cgal_program( "${cppfile}" )
 # endforeach()
 #todo later: re-add this foreach (after the merge with master, whith the last viewers)
 create_single_source_cgal_program(map_2_constructor.cpp)

 add_executable(path_homotopy path_homotopy.cpp)
 target_link_libraries(path_homotopy ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})

  if (USE_BASIC_VIEWER)
    target_link_libraries(path_homotopy ${BASIC_VIEWER_LIBRARIES})
    qt5_use_modules(path_homotopy ${BASIC_VIEWER_MODULES})
  endif(USE_BASIC_VIEWER)

else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()
