# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.12...3.31)
project(Surface_mesh_deformation_Examples)

find_package(CGAL REQUIRED)

find_package(Eigen3 3.1.91 QUIET) #(requires 3.1.91 or greater)
include(CGAL_Eigen3_support)
if(TARGET CGAL::Eigen3_support)
  create_single_source_cgal_program("all_roi_assign_example.cpp")
  create_single_source_cgal_program("all_roi_assign_example_custom_polyhedron.cpp")
  create_single_source_cgal_program("all_roi_assign_example_Surface_mesh.cpp")
  create_single_source_cgal_program("custom_weight_for_edges_example.cpp")
  create_single_source_cgal_program("deform_polyhedron_with_custom_pmap_example.cpp")
  create_single_source_cgal_program("k_ring_roi_translate_rotate_example.cpp")
  create_single_source_cgal_program("k_ring_roi_translate_rotate_Surface_mesh.cpp")
  create_single_source_cgal_program("deform_mesh_for_botsch08_format_sre_arap.cpp")

  foreach(
    target
    all_roi_assign_example
    all_roi_assign_example_custom_polyhedron
    all_roi_assign_example_Surface_mesh
    custom_weight_for_edges_example
    deform_polyhedron_with_custom_pmap_example
    k_ring_roi_translate_rotate_example
    k_ring_roi_translate_rotate_Surface_mesh
    deform_mesh_for_botsch08_format_sre_arap)
    target_link_libraries(${target} PRIVATE CGAL::Eigen3_support)
  endforeach()

  find_package(OpenMesh QUIET)
  if(OpenMesh_FOUND)
    include(CGAL_OpenMesh_support)
    message(STATUS "Found OpenMesh")

    create_single_source_cgal_program("all_roi_assign_example_with_OpenMesh.cpp")
    target_link_libraries(all_roi_assign_example_with_OpenMesh
                          PRIVATE CGAL::OpenMesh_support CGAL::Eigen3_support)
  else()
    message(STATUS "NOTICE: Examples that use OpenMesh will not be compiled.")
  endif()

else()
  message("NOTICE: These examples require Eigen 3.1.91 (or greater) or later and will not be compiled.")
endif()
