# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.12...3.29)
project(BGL_polyhedron_3_Examples)

# CGAL and its components
find_package(CGAL REQUIRED)

create_single_source_cgal_program("distance.cpp")
create_single_source_cgal_program("incident_vertices.cpp")
create_single_source_cgal_program("kruskal.cpp")
create_single_source_cgal_program("kruskal_with_stored_id.cpp")
create_single_source_cgal_program("normals.cpp")
create_single_source_cgal_program("range.cpp")
create_single_source_cgal_program("transform_iterator.cpp")
create_single_source_cgal_program("copy_polyhedron.cpp")

find_package(OpenMesh QUIET)
if(OpenMesh_FOUND)
  include(CGAL_OpenMesh_support)
  target_link_libraries(copy_polyhedron PRIVATE CGAL::OpenMesh_support)
else()
  message(STATUS "NOTICE: The example 'copy_polyhedron' will not use OpenMesh.")
endif()

find_package(METIS QUIET)
include(CGAL_METIS_support)
if(TARGET CGAL::METIS_support)
  create_single_source_cgal_program("polyhedron_partition.cpp")
  target_link_libraries(polyhedron_partition PUBLIC CGAL::METIS_support)
else()
  message(STATUS "NOTICE: The example 'polyhedron_partition' requires the METIS library, and will not be compiled.")
endif()
