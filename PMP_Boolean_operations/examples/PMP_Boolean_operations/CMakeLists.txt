# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.12...3.31)
project(PMP_Boolean_operations_Examples)

# CGAL and its components
find_package(CGAL REQUIRED)

create_single_source_cgal_program("corefinement_SM.cpp")
create_single_source_cgal_program("corefinement_consecutive_bool_op.cpp")
create_single_source_cgal_program("corefinement_difference_remeshed.cpp")
create_single_source_cgal_program("corefinement_mesh_union.cpp")
create_single_source_cgal_program("corefinement_mesh_non_manifold_intersection.cpp")
create_single_source_cgal_program("corefinement_mesh_union_progress.cpp")
create_single_source_cgal_program("corefinement_mesh_union_and_intersection.cpp")
create_single_source_cgal_program("corefinement_mesh_union_with_attributes.cpp")
create_single_source_cgal_program("corefinement_polyhedron_union.cpp")
create_single_source_cgal_program("corefinement_LCC.cpp")
create_single_source_cgal_program("soup_autorefinement.cpp")
create_single_source_cgal_program("surface_mesh_intersection.cpp")
create_single_source_cgal_program("triangle_mesh_autorefinement.cpp")
create_single_source_cgal_program("snap_polygon_soup.cpp")
create_single_source_cgal_program("mesh_slicer_example.cpp")
create_single_source_cgal_program("mesh_kernel_example.cpp")

find_package(Eigen3 QUIET)
include(CGAL_Eigen3_support)

if(TARGET CGAL::Eigen3_support)
  create_single_source_cgal_program("geodesic_isolevel_refinement.cpp")
  target_link_libraries(geodesic_isolevel_refinement PRIVATE CGAL::Eigen3_support)
else()
  message(STATUS "NOTICE: Examples that use Eigen will not be compiled.")
endif()

find_package(OpenMesh QUIET)
if(OpenMesh_FOUND)
  include(CGAL_OpenMesh_support)

  create_single_source_cgal_program("corefinement_OM_union.cpp")
  target_link_libraries(corefinement_OM_union PRIVATE CGAL::OpenMesh_support)
else()
  message(STATUS "NOTICE: Examples that use OpenMesh will not be compiled.")
endif()

find_package(TBB QUIET)
include(CGAL_TBB_support)
if(TARGET CGAL::TBB_support)
  target_link_libraries(soup_autorefinement PRIVATE CGAL::TBB_support)
  target_link_libraries(snap_polygon_soup PRIVATE CGAL::TBB_support)

  create_single_source_cgal_program("corefinement_parallel_union_meshes.cpp")
  target_link_libraries(corefinement_parallel_union_meshes PRIVATE CGAL::TBB_support)
else()
  message(STATUS "NOTICE: Intel TBB was not found. Sequential code will be used.")
endif()

