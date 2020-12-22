# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.15)
project(Surface_mesh_deformation_Demo)

if(NOT POLICY CMP0070 AND POLICY CMP0053)
  # Only set CMP0053 to OLD with CMake<3.10, otherwise there is a warning.
  cmake_policy(SET CMP0053 OLD)
endif()

set_property(DIRECTORY PROPERTY CGAL_NO_TESTING TRUE)

find_package(CGAL REQUIRED)

find_package(Eigen3 3.1.91) #(requires 3.2.0 or greater)
include(CGAL_Eigen_support)
if(TARGET CGAL::Eigen_support)
  create_single_source_cgal_program("deform_mesh_for_botsch08_format.cpp")
  target_link_libraries(deform_mesh_for_botsch08_format
                        PUBLIC CGAL::Eigen_support)
else()
  message(
    STATUS
      "NOTICE: This program requires the Eigen library, version 3.2 or later and will not be compiled."
  )
endif()
