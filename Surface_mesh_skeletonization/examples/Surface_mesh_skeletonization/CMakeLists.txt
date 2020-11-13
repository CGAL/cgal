# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.15)
project(Surface_mesh_skeletonization_Examples)

find_package(CGAL REQUIRED)

find_package(Eigen3 3.2.0) #(requires 3.2.0 or greater)
include(CGAL_Eigen_support)

if(TARGET CGAL::Eigen_support)
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
    target_link_libraries(${target} PUBLIC CGAL::Eigen_support)
  endforeach()
else()
  message(
    STATUS
      "These programs require the Eigen library (3.2 or greater), and will not be compiled."
  )
endif()
