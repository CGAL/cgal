# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.1...3.15)
project(Heat_method_3_Examples)

# CGAL and its components
find_package(CGAL REQUIRED)

# include helper file
include(${CGAL_USE_FILE})

# Boost and its components
find_package(Boost REQUIRED)

if(NOT Boost_FOUND)

  message(
    STATUS "This project requires the Boost library, and will not be compiled.")

  return()

endif()

find_package(Eigen3 3.3.0)
include(CGAL_Eigen_support)

if(NOT TARGET CGAL::Eigen_support)
  message(
    STATUS
      "This project requires the Eigen library (3.3 or greater), and will not be compiled."
  )
  return()
endif()

# include for local directory
include_directories(BEFORE include)

# Creating entries for all C++ files with "main" routine
# ##########################################################

include(CGAL_CreateSingleSourceCGALProgram)

create_single_source_cgal_program("heat_method.cpp")
target_link_libraries(heat_method PUBLIC CGAL::Eigen_support)
create_single_source_cgal_program("heat_method_polyhedron.cpp")
target_link_libraries(heat_method_polyhedron PUBLIC CGAL::Eigen_support)
create_single_source_cgal_program("heat_method_surface_mesh.cpp")
target_link_libraries(heat_method_surface_mesh PUBLIC CGAL::Eigen_support)
create_single_source_cgal_program("heat_method_surface_mesh_direct.cpp")
target_link_libraries(heat_method_surface_mesh_direct
                      PUBLIC CGAL::Eigen_support)
