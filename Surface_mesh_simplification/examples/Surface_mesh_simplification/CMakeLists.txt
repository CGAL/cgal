# Created by the script cgal_create_cmake_script_with_options
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.1...3.15)
project(Surface_mesh_simplification_Examples)

# CGAL and its components
find_package(CGAL REQUIRED)

include(${CGAL_USE_FILE})

# Boost and its components
find_package(Boost)
if(NOT Boost_FOUND)
  message(
    STATUS "This project requires the Boost library, and will not be compiled.")
  return()
endif()

find_package(OpenMesh QUIET)

if(OpenMesh_FOUND)
  include(UseOpenMesh)
else()
  message(STATUS "Examples that use OpenMesh will not be compiled.")
endif()

# Creating entries for all .cpp/.C files with "main" routine
# ##########################################################

create_single_source_cgal_program("edge_collapse_constrain_sharp_edges.cpp")
create_single_source_cgal_program(
  "edge_collapse_constrained_border_polyhedron.cpp")
create_single_source_cgal_program("edge_collapse_enriched_polyhedron.cpp")
create_single_source_cgal_program("edge_collapse_polyhedron.cpp")
create_single_source_cgal_program("edge_collapse_surface_mesh.cpp")
create_single_source_cgal_program("edge_collapse_linear_cell_complex.cpp")
create_single_source_cgal_program(
  "edge_collapse_constrained_border_surface_mesh.cpp")
create_single_source_cgal_program("edge_collapse_all_short_edges.cpp")
create_single_source_cgal_program("edge_collapse_bounded_normal_change.cpp")
create_single_source_cgal_program("edge_collapse_visitor_surface_mesh.cpp")

find_package(Eigen3 3.1.0 QUIET) #(3.1.0 or greater)
include(CGAL_Eigen_support)
if(TARGET CGAL::Eigen_support)
  create_single_source_cgal_program("edge_collapse_garland_heckbert.cpp")
  target_link_libraries(edge_collapse_garland_heckbert
                        PUBLIC CGAL::Eigen_support)
else()
  message(
    STATUS
      "Garland-Heckbert polices require the Eigen library, which has not been found; related examples will not be compiled."
  )
endif()

if(OpenMesh_FOUND)
  create_single_source_cgal_program("edge_collapse_OpenMesh.cpp")
  target_link_libraries(edge_collapse_OpenMesh PRIVATE ${OPENMESH_LIBRARIES})
endif()
