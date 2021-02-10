# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.15)
project(Surface_mesh_shortest_path_Examples)

find_package(CGAL REQUIRED)

create_single_source_cgal_program("shortest_path_sequence.cpp")
create_single_source_cgal_program("shortest_paths_multiple_sources.cpp")
create_single_source_cgal_program("shortest_paths_no_id.cpp")
create_single_source_cgal_program("shortest_paths_with_id.cpp")
create_single_source_cgal_program("shortest_paths.cpp")
create_single_source_cgal_program("shortest_path_with_locate.cpp")

find_package(OpenMesh QUIET)
if(OpenMesh_FOUND)
  include(UseOpenMesh)
  create_single_source_cgal_program("shortest_paths_OpenMesh.cpp")
  target_link_libraries(shortest_paths_OpenMesh PRIVATE ${OPENMESH_LIBRARIES})
else()
  message(STATUS "Examples that use OpenMesh will not be compiled.")
endif()
