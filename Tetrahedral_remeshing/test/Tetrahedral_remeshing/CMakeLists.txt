# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.1...3.23)

project(Tetrahedral_remeshing_Tests)

# CGAL and its components
find_package(CGAL REQUIRED COMPONENTS ImageIO)

create_single_source_cgal_program("test_tetrahedral_remeshing.cpp")
create_single_source_cgal_program("test_tetrahedral_remeshing_with_features.cpp")
create_single_source_cgal_program("test_tetrahedral_remeshing_of_one_subdomain.cpp")
create_single_source_cgal_program("test_tetrahedral_remeshing_io.cpp")
create_single_source_cgal_program("test_tetrahedral_remeshing_from_mesh_file.cpp")

# Tests using Mesh_3 require Eigen
find_package(Eigen3 3.1.0 QUIET) #(3.1.0 or greater)
include(CGAL_Eigen3_support)
if(TARGET CGAL::Eigen3_support)
  create_single_source_cgal_program("test_mesh_and_remesh_polyhedral_domain_with_features.cpp")
  target_link_libraries(test_mesh_and_remesh_polyhedral_domain_with_features PUBLIC CGAL::Eigen3_support)

  if(CGAL_ImageIO_USE_ZLIB)
    create_single_source_cgal_program("test_mesh_and_remesh_image.cpp")
    target_link_libraries(test_mesh_and_remesh_image PUBLIC CGAL::Eigen3_support)
  else()
    message(STATUS "NOTICE: The test 'test_mesh_and_remesh_image' uses zlib, and will not be compiled.")
  endif()
else()
  message(STATUS "NOTICE: Some tests require Eigen 3.1 (or greater), and will not be compiled.")
endif()
