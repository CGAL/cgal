# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.12...3.31)
project(Surface_mesh_skeletonization_Examples)

find_package(CGAL REQUIRED)

find_package(Eigen3 QUIET)
include(CGAL_Eigen3_support)

if(TARGET CGAL::Eigen3_support)
  create_single_source_cgal_program("simple_mcfskel_example.cpp")
  create_single_source_cgal_program("simple_mcfskel_sm_example.cpp")
  create_single_source_cgal_program("simple_mcfskel_LCC_example.cpp")
  create_single_source_cgal_program("MCF_Skeleton_example.cpp")
  create_single_source_cgal_program("MCF_Skeleton_sm_example.cpp")
  create_single_source_cgal_program("MCF_Skeleton_LCC_example.cpp")
  create_single_source_cgal_program("segmentation_example.cpp")

  foreach(
    target
    simple_mcfskel_example
    simple_mcfskel_sm_example
    simple_mcfskel_LCC_example
    MCF_Skeleton_example
    MCF_Skeleton_sm_example
    MCF_Skeleton_LCC_example
    segmentation_example)
    target_link_libraries(${target} PRIVATE CGAL::Eigen3_support)
  endforeach()

  find_package(OpenMesh QUIET)
  if(OpenMesh_FOUND)
    include(CGAL_OpenMesh_support)
    message(STATUS "Found OpenMesh")

    create_single_source_cgal_program("MCF_Skeleton_om_example.cpp")
    target_link_libraries(MCF_Skeleton_om_example PRIVATE CGAL::Eigen3_support CGAL::OpenMesh_support)
  else()
    message(STATUS "NOTICE: Examples that use OpenMesh will not be compiled.")
  endif()
else()
  message("NOTICE: These programs require the Eigen library (3.2 or greater), and will not be compiled.")
endif()
