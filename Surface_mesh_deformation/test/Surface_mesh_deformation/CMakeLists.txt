# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.20)
project(Surface_mesh_deformation_Tests)

find_package(CGAL REQUIRED)

find_package(Eigen3 3.1.91) #(requires 3.2.0 or greater)
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
    include(UseOpenMesh)
    create_single_source_cgal_program("Cactus_deformation_session_OpenMesh.cpp")
    target_link_libraries(Cactus_deformation_session_OpenMesh
                          PRIVATE ${OPENMESH_LIBRARIES} CGAL::Eigen3_support)
  else()
    message(STATUS "Example that use OpenMesh will not be compiled.")
  endif()
else()
  message(
    STATUS
      "NOTICE: These tests require the Eigen library, version 3.2 or later and will not be compiled."
  )
endif()
