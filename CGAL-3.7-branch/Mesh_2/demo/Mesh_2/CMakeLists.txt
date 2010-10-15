# Created by the script cgal_create_cmake_script (and then adapted manually).
# This is the CMake script for compiling a CGAL application.

project( mesh_2_demo )

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

add_subdirectory(Qt3)
 
if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )  
endif()

find_package(CGAL COMPONENTS Core )
include( ${CGAL_USE_FILE} )

include_directories(BEFORE ../../include)
include_directories(BEFORE ../../../Triangulation_2/include)


if ( CGAL_FOUND )

  # conform target
  add_executable  (conform conform.cpp)
  target_link_libraries(conform ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )
  add_to_cached_list( CGAL_EXECUTABLE_TARGETS conform )

  # mesh target
  add_executable  (mesh mesh.cpp)
  target_link_libraries(mesh ${QT3_LIBRARIES} ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )
  add_to_cached_list( CGAL_EXECUTABLE_TARGETS mesh )

else()

  message(STATUS "NOTICE: This demo requires CGAL and will not be compiled.")

endif()
