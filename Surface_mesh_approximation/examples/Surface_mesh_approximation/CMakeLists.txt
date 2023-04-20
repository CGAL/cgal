# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.1...3.23)
project(Surface_mesh_approximation_Examples)

# CGAL and its components
find_package(CGAL REQUIRED)

# Use Eigen (for PCA)
find_package(Eigen3 3.1.0 QUIET) #(3.1.0 or greater)
include(CGAL_Eigen3_support)
if(NOT TARGET CGAL::Eigen3_support)
  message("NOTICE: This project requires Eigen 3.1 (or greater), and will not be compiled.")
  return()
endif()

create_single_source_cgal_program("vsa_approximation_2_example.cpp")
target_link_libraries(vsa_approximation_2_example PUBLIC CGAL::Eigen3_support)

create_single_source_cgal_program("vsa_approximation_example.cpp")
target_link_libraries(vsa_approximation_example PUBLIC CGAL::Eigen3_support)

create_single_source_cgal_program("vsa_class_interface_example.cpp")
target_link_libraries(vsa_class_interface_example PUBLIC CGAL::Eigen3_support)

create_single_source_cgal_program("vsa_isotropic_metric_example.cpp")
target_link_libraries(vsa_isotropic_metric_example PUBLIC CGAL::Eigen3_support)

create_single_source_cgal_program("vsa_segmentation_example.cpp")
target_link_libraries(vsa_segmentation_example PUBLIC CGAL::Eigen3_support)

create_single_source_cgal_program("vsa_simple_approximation_example.cpp")
target_link_libraries(vsa_simple_approximation_example PUBLIC CGAL::Eigen3_support)
