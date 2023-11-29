# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.1...3.23)
project(Surface_mesh_approximation_Tests)

# CGAL and its components
find_package(CGAL REQUIRED)

# Use Eigen (for PCA)
find_package(Eigen3 3.1.0 QUIET) #(3.1.0 or greater)
include(CGAL_Eigen3_support)
if(NOT TARGET CGAL::Eigen3_support)
  message("NOTICE: This project requires Eigen 3.1 (or greater), and will not be compiled.")
  return()
endif()

create_single_source_cgal_program("vsa_border_test.cpp")
target_link_libraries(vsa_border_test PUBLIC CGAL::Eigen3_support)

create_single_source_cgal_program("vsa_class_interface_test.cpp")
target_link_libraries(vsa_class_interface_test PUBLIC CGAL::Eigen3_support)

create_single_source_cgal_program("vsa_correctness_test.cpp")
target_link_libraries(vsa_correctness_test PUBLIC CGAL::Eigen3_support)

create_single_source_cgal_program("vsa_error_decrease_test.cpp")
target_link_libraries(vsa_error_decrease_test PUBLIC CGAL::Eigen3_support)

create_single_source_cgal_program("vsa_kernel_test.cpp")
target_link_libraries(vsa_kernel_test PUBLIC CGAL::Eigen3_support)

create_single_source_cgal_program("vsa_approximation_test.cpp")
target_link_libraries(vsa_approximation_test PUBLIC CGAL::Eigen3_support)

create_single_source_cgal_program("vsa_segmentation_test.cpp")
target_link_libraries(vsa_segmentation_test PUBLIC CGAL::Eigen3_support)

#create_single_source_cgal_program("vsa_meshing_manifold_test.cpp")
#target_link_libraries(vsa_meshing_manifold_test PUBLIC CGAL::Eigen3_support)

create_single_source_cgal_program("vsa_metric_test.cpp")
target_link_libraries(vsa_metric_test PUBLIC CGAL::Eigen3_support)

create_single_source_cgal_program("vsa_teleportation_test.cpp")
target_link_libraries(vsa_teleportation_test PUBLIC CGAL::Eigen3_support)
