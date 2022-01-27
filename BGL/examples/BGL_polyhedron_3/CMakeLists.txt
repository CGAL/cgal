# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.1...3.22)
project(BGL_polyhedron_3_Examples)

# CGAL and its components
find_package(CGAL REQUIRED)

# Boost and its components
find_package(Boost REQUIRED)

if(NOT Boost_FOUND)

  message(
    STATUS "This project requires the Boost library, and will not be compiled.")

  return()

endif()

# include for local directory

# include for local package

# Creating entries for all C++ files with "main" routine
# ##########################################################

create_single_source_cgal_program("distance.cpp")

create_single_source_cgal_program("incident_vertices.cpp")

create_single_source_cgal_program("kruskal.cpp")

create_single_source_cgal_program("kruskal_with_stored_id.cpp")

create_single_source_cgal_program("normals.cpp")

create_single_source_cgal_program("range.cpp")

create_single_source_cgal_program("transform_iterator.cpp")

create_single_source_cgal_program("copy_polyhedron.cpp")

find_package( OpenMesh QUIET )
if(OpenMesh_FOUND)
  target_link_libraries( copy_polyhedron PRIVATE ${OPENMESH_LIBRARIES} )
  target_compile_definitions( copy_polyhedron PRIVATE -DCGAL_USE_OPENMESH )
else()
  message(STATUS "Examples that use OpenMesh will not be compiled.")
endif()

find_package( METIS )
include(CGAL_METIS_support)
if( TARGET CGAL::METIS_support )
  create_single_source_cgal_program( "polyhedron_partition.cpp" )
  target_link_libraries( polyhedron_partition PUBLIC CGAL::METIS_support)
else()
  message( STATUS "Examples that use the METIS library will not be compiled." )
endif()
