# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.12...3.31)
project(PMP_Boolean_operations_Tests)

# CGAL and its components
find_package(CGAL REQUIRED COMPONENTS Core)

create_single_source_cgal_program("autorefinement_sm.cpp")
create_single_source_cgal_program("corefine_non_manifold.cpp" )
create_single_source_cgal_program("polygon_mesh_slicer_test.cpp")
create_single_source_cgal_program("surface_intersection_sm_poly.cpp")
create_single_source_cgal_program("test_autorefinement.cpp")
create_single_source_cgal_program("test_corefine.cpp")
create_single_source_cgal_program("test_coref_epic_points_identity.cpp")
create_single_source_cgal_program("test_corefinement_and_constraints.cpp")
create_single_source_cgal_program("test_corefinement_and_constraints_nearest.cpp")
create_single_source_cgal_program("test_corefinement_bool_op.cpp")
create_single_source_cgal_program("test_corefinement_cavities.cpp")
create_single_source_cgal_program("test_corefinement_nm_bo.cpp")
create_single_source_cgal_program("test_degenerate_pmp_clip_split_corefine.cpp")
create_single_source_cgal_program("test_isolevel_refinement.cpp")
create_single_source_cgal_program("test_pmp_clip.cpp")
create_single_source_cgal_program("test_snap_rounding.cpp")
create_single_source_cgal_program("test_mesh_kernel.cpp")

find_package(TBB QUIET)
include(CGAL_TBB_support)
if(TARGET CGAL::TBB_support)
  target_link_libraries(test_autorefinement PRIVATE CGAL::TBB_support)
  target_link_libraries(test_snap_rounding PRIVATE CGAL::TBB_support)
else()
  message(STATUS "NOTICE: Intel TBB was not found. Tests will use sequential code.")
endif()
