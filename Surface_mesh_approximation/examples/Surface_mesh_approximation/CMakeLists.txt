# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.1...3.15)
project(Surface_mesh_approximation_Examples)

# CGAL and its components
find_package(CGAL REQUIRED)

# include helper file
include(${CGAL_USE_FILE})

# Boost
find_package(Boost)
if(NOT Boost_FOUND)
  message(
    STATUS "This project requires the Boost library, and will not be compiled.")
  return()
endif()

# Use Eigen (for PCA)
find_package(Eigen3 3.1.0 REQUIRED) #(3.1.0 or greater)
include(CGAL_Eigen_support)
if(NOT TARGET CGAL::Eigen_support)
  message(
    STATUS "This project requires the Eigen library, and will not be compiled.")
  return()
endif()

# Creating entries for all C++ files with "main" routine
# ##########################################################

include(CGAL_CreateSingleSourceCGALProgram)

create_single_source_cgal_program("vsa_approximation_2_example.cpp")
target_link_libraries(vsa_approximation_2_example PUBLIC CGAL::Eigen_support)

create_single_source_cgal_program("vsa_approximation_example.cpp")
target_link_libraries(vsa_approximation_example PUBLIC CGAL::Eigen_support)

create_single_source_cgal_program("vsa_class_interface_example.cpp")
target_link_libraries(vsa_class_interface_example PUBLIC CGAL::Eigen_support)

create_single_source_cgal_program("vsa_isotropic_metric_example.cpp")
target_link_libraries(vsa_isotropic_metric_example PUBLIC CGAL::Eigen_support)

create_single_source_cgal_program("vsa_segmentation_example.cpp")
target_link_libraries(vsa_segmentation_example PUBLIC CGAL::Eigen_support)

create_single_source_cgal_program("vsa_simple_approximation_example.cpp")
target_link_libraries(vsa_simple_approximation_example
                      PUBLIC CGAL::Eigen_support)
