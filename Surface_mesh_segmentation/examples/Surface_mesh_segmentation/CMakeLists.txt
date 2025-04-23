# Created by the script cgal_create_cmake_script_with_options
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.12...3.31)
project(Surface_mesh_segmentation_Examples)

# CGAL and its components
find_package(CGAL REQUIRED)

create_single_source_cgal_program("sdf_values_example.cpp")
create_single_source_cgal_program("segmentation_from_sdf_values_example.cpp")
create_single_source_cgal_program("segmentation_via_sdf_values_example.cpp")
create_single_source_cgal_program("segmentation_with_facet_ids_example.cpp")
create_single_source_cgal_program("segmentation_from_sdf_values_SM_example.cpp")
create_single_source_cgal_program("segmentation_from_sdf_values_LCC_example.cpp")
create_single_source_cgal_program("extract_segmentation_into_mesh_example.cpp")

find_package(OpenMesh QUIET)
if(OpenMesh_FOUND)
  include(CGAL_OpenMesh_support)
  message(STATUS "Found OpenMesh")

  create_single_source_cgal_program("segmentation_from_sdf_values_OpenMesh_example.cpp")
  target_link_libraries(segmentation_from_sdf_values_OpenMesh_example PRIVATE CGAL::OpenMesh_support)
else()
  message(STATUS "NOTICE: Examples that use OpenMesh will not be compiled.")
endif()
