# Created by the script cgal_create_cmake_script (and then adapted manually).
# This is the CMake script for compiling a CGAL application.

project( mesh_2_demo )

cmake_minimum_required(VERSION 2.8.11)
if(POLICY CMP0043)
  cmake_policy(SET CMP0043 OLD)
endif()

add_subdirectory(Qt3)
 
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
  target_link_libraries(mesh ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )
  add_to_cached_list( CGAL_EXECUTABLE_TARGETS mesh )

else()

  message(STATUS "NOTICE: This demo requires CGAL and will not be compiled.")

endif()
