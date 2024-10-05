# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.12...3.29)
project(Surface_mesh_deformation_Tests)

find_package(CGAL REQUIRED)

find_package(Eigen3 3.1.91 QUIET) #(requires 3.1.91 or greater)
include(CGAL_Eigen3_support)
if(TARGET CGAL::Eigen3_support)
  create_single_source_cgal_program("Cactus_deformation_session.cpp")
  target_link_libraries(Cactus_deformation_session PUBLIC CGAL::Eigen3_support)
  create_single_source_cgal_program("Cactus_performance_test.cpp")
  target_link_libraries(Cactus_performance_test PUBLIC CGAL::Eigen3_support)
  create_single_source_cgal_program("Symmetry_test.cpp")
  target_link_libraries(Symmetry_test PUBLIC CGAL::Eigen3_support)

  find_package(OpenMesh QUIET)
  if(OpenMesh_FOUND)
    include(CGAL_OpenMesh_support)
    message(STATUS "Found OpenMesh")

    create_single_source_cgal_program("Cactus_deformation_session_OpenMesh.cpp")
    target_link_libraries(Cactus_deformation_session_OpenMesh
                          PRIVATE CGAL::OpenMesh_support CGAL::Eigen3_support)
  else()
    message(STATUS "NOTICE: Examples that use OpenMesh will not be compiled.")
  endif()
else()
  message("NOTICE: These tests require the Eigen library, version 3.1.91 or later and will not be compiled.")
endif()
