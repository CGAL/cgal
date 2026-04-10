# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.12...3.31)
project(PMP_Remeshing_Tests)

# CGAL and its components
find_package(CGAL REQUIRED COMPONENTS Core)

create_single_source_cgal_program("extrude_test.cpp")
create_single_source_cgal_program("remeshing_test.cpp")
create_single_source_cgal_program("remeshing_test_epeck.cpp")
create_single_source_cgal_program("remeshing_with_isolated_constraints_test.cpp" )
create_single_source_cgal_program("test_mesh_smoothing.cpp")
create_single_source_cgal_program("triangulate_faces_test.cpp")
create_single_source_cgal_program("triangulate_faces_hole_filling_dt3_test.cpp")
create_single_source_cgal_program("triangulate_faces_hole_filling_all_search_test.cpp")

find_package(Eigen3 QUIET)
include(CGAL_Eigen3_support)

if(TARGET CGAL::Eigen3_support)
  create_single_source_cgal_program("delaunay_remeshing_test.cpp")
  target_link_libraries(delaunay_remeshing_test PRIVATE CGAL::Eigen3_support)
  create_single_source_cgal_program("fairing_test.cpp")
  target_link_libraries(fairing_test PRIVATE CGAL::Eigen3_support)
  create_single_source_cgal_program("test_shape_smoothing.cpp")
  target_link_libraries(test_shape_smoothing PRIVATE CGAL::Eigen3_support)
  create_single_source_cgal_program("test_decimation_of_planar_patches.cpp")
  target_link_libraries(test_decimation_of_planar_patches PRIVATE CGAL::Eigen3_support)
  create_single_source_cgal_program("remeshing_quality_test.cpp" )
  target_link_libraries(remeshing_quality_test PRIVATE CGAL::Eigen3_support)
  create_single_source_cgal_program("test_acvd.cpp" )
  target_link_libraries(test_acvd PRIVATE CGAL::Eigen3_support)
else()
  message(STATUS "NOTICE: Tests that use the Eigen library will not be compiled.")
endif()

find_package(OpenMesh QUIET)
if(OpenMesh_FOUND)
  include(CGAL_OpenMesh_support)
  create_single_source_cgal_program("remeshing_test_P_SM_OM.cpp")
  target_link_libraries(remeshing_test_P_SM_OM PRIVATE CGAL::OpenMesh_support)
else()
  message(STATUS "NOTICE: Tests that use OpenMesh will not be compiled.")
endif()

find_package(Ceres QUIET)
include(CGAL_Ceres_support)
if(TARGET CGAL::Ceres_support AND TARGET CGAL::Eigen3_support)
  target_link_libraries(test_mesh_smoothing PRIVATE CGAL::Eigen3_support CGAL::Ceres_support)
else()
  message(STATUS "NOTICE: Tests are not using Ceres.")
endif()
