# Created by the script cgal_create_cmake_script.
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.12...3.29)
project(Shape_detection_Tests)

find_package(CGAL REQUIRED COMPONENTS Core)

create_single_source_cgal_program("test_efficient_RANSAC_cone_connected_component.cpp")
create_single_source_cgal_program("test_efficient_RANSAC_cone_parameters.cpp")
create_single_source_cgal_program("test_efficient_RANSAC_cylinder_connected_component.cpp")
create_single_source_cgal_program("test_efficient_RANSAC_cylinder_parameters.cpp")
create_single_source_cgal_program("test_efficient_RANSAC_plane_connected_component.cpp")
create_single_source_cgal_program("test_efficient_RANSAC_plane_parameters.cpp")
create_single_source_cgal_program("test_efficient_RANSAC_sphere_connected_component.cpp")
create_single_source_cgal_program("test_efficient_RANSAC_sphere_parameters.cpp")
create_single_source_cgal_program("test_efficient_RANSAC_torus_connected_component.cpp")
create_single_source_cgal_program("test_efficient_RANSAC_torus_parameters.cpp")
create_single_source_cgal_program("test_efficient_RANSAC_scene.cpp")

find_package(Eigen3 3.1.0 QUIET) # (3.1.0 or greater)
include(CGAL_Eigen3_support)
if(TARGET CGAL::Eigen3_support)

  create_single_source_cgal_program("test_region_growing_basic.cpp")
  create_single_source_cgal_program("test_region_growing_strict.cpp")
  create_single_source_cgal_program("test_region_growing_on_cube.cpp")
  create_single_source_cgal_program("test_region_growing_on_point_set_2.cpp")
  create_single_source_cgal_program("test_region_growing_on_point_set_3.cpp")
  create_single_source_cgal_program("test_region_growing_on_polygon_mesh.cpp")
  create_single_source_cgal_program("test_region_growing_on_point_set_2_with_sorting.cpp")
  create_single_source_cgal_program("test_region_growing_on_point_set_3_with_sorting.cpp")
  create_single_source_cgal_program("test_region_growing_on_polygon_mesh_with_sorting.cpp")
  create_single_source_cgal_program("test_region_growing_on_degenerated_mesh.cpp")

  foreach(
    target
    test_region_growing_basic
    test_region_growing_strict
    test_region_growing_on_cube
    test_region_growing_on_point_set_2
    test_region_growing_on_point_set_3
    test_region_growing_on_polygon_mesh
    test_region_growing_on_point_set_2_with_sorting
    test_region_growing_on_point_set_3_with_sorting
    test_region_growing_on_polygon_mesh_with_sorting
    test_region_growing_on_degenerated_mesh)
    target_link_libraries(${target} PUBLIC CGAL::Eigen3_support)
  endforeach()

  set(RANSAC_PROTO_DIR CACHE PATH "")
  if(NOT RANSAC_PROTO_DIR STREQUAL "")
    add_definitions(-DPOINTSWITHINDEX -DCGAL_TEST_RANSAC_PROTOTYPE)
    include_directories(${RANSAC_PROTO_DIR})
    include_directories(${RANSAC_PROTO_DIR}/MiscLib/)
    file(GLOB proto_src "${RANSAC_PROTO_DIR}/*.cpp")
    file(GLOB proto_misc_src "${RANSAC_PROTO_DIR}/MiscLib/*.cpp")
    add_library(libproto STATIC  ${proto_src} ${proto_misc_src})
    add_executable(test_validity_sampled_data "test_validity_sampled_data.cpp")
    target_link_libraries(test_validity_sampled_data libproto CGAL::CGAL CGAL::Data CGAL::Eigen3_support)
  else()
    add_executable(test_validity_sampled_data "test_validity_sampled_data.cpp")
    target_link_libraries(test_validity_sampled_data CGAL::CGAL CGAL::Data CGAL::Eigen3_support)
  endif()
  cgal_add_test(test_validity_sampled_data)
else()
  message(STATUS "NOTICE: Some tests require Eigen 3.1 (or greater), and will not be compiled.")
endif()
