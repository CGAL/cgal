# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.12...3.31)

project(Tetrahedral_remeshing_Examples)

# CGAL and its components
find_package(CGAL REQUIRED)

create_single_source_cgal_program("tetrahedral_remeshing_example.cpp" )
create_single_source_cgal_program("tetrahedral_remeshing_with_features.cpp")
create_single_source_cgal_program("tetrahedral_remeshing_of_one_subdomain.cpp")
create_single_source_cgal_program("tetrahedral_remeshing_from_mesh.cpp")

# Concurrent Mesh_3
option(CGAL_ACTIVATE_CONCURRENT_MESH_3 "Activate parallelism in Mesh_3" OFF)
if(CGAL_ACTIVATE_CONCURRENT_MESH_3 OR "$ENV{CGAL_ACTIVATE_CONCURRENT_MESH_3}")
  add_compile_definitions(CGAL_CONCURRENT_MESH_3)
  find_package(TBB REQUIRED)
  include(CGAL_TBB_support)
else()
  find_package(TBB QUIET)
  include(CGAL_TBB_support)
endif()

# Use Eigen for Mesh_3
find_package(Eigen3 QUIET)
include(CGAL_Eigen3_support)

if(TARGET CGAL::Eigen3_support)
  create_single_source_cgal_program( "mesh_and_remesh_polyhedral_domain_with_features.cpp" )
  target_link_libraries(mesh_and_remesh_polyhedral_domain_with_features PRIVATE CGAL::Eigen3_support)

  create_single_source_cgal_program("mesh_and_remesh_c3t3.cpp")
  target_link_libraries(mesh_and_remesh_c3t3 PRIVATE CGAL::Eigen3_support)

  create_single_source_cgal_program( "mesh_and_remesh_with_adaptive_sizing.cpp")
  target_link_libraries(mesh_and_remesh_with_adaptive_sizing PRIVATE CGAL::Eigen3_support)

  create_single_source_cgal_program( "mesh_and_remesh_with_sizing.cpp")
  target_link_libraries(mesh_and_remesh_with_sizing PRIVATE CGAL::Eigen3_support)

  create_single_source_cgal_program("tetrahedral_remeshing_with_features_and_sizing.cpp")
  target_link_libraries(tetrahedral_remeshing_with_features_and_sizing PRIVATE CGAL::Eigen3_support)

  if(CGAL_ACTIVATE_CONCURRENT_MESH_3 AND TARGET CGAL::TBB_support)
    message(STATUS "Found TBB")
    target_link_libraries(mesh_and_remesh_polyhedral_domain_with_features PRIVATE CGAL::TBB_support)
    target_link_libraries(mesh_and_remesh_c3t3 PRIVATE CGAL::TBB_support)
    target_link_libraries(mesh_and_remesh_with_sizing PRIVATE CGAL::TBB_support)
    target_link_libraries(tetrahedral_remeshing_with_features_and_sizing PRIVATE CGAL::TBB_support)
  endif()
else()
  message(STATUS "NOTICE: Some examples require Eigen 3.1 (or greater), and will not be compiled.")
endif()

