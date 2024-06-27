# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.12...3.29)

project(Tetrahedral_remeshing_Tests)

# CGAL and its components
find_package(CGAL REQUIRED COMPONENTS ImageIO)

add_definitions(-DCGAL_TETRAHEDRAL_REMESHING_VERBOSE)
#add_definitions(-DCGAL_CHECK_EXPENSIVE) #to activate expensive assertions

create_single_source_cgal_program("test_tetrahedral_remeshing.cpp")
create_single_source_cgal_program("test_tetrahedral_remeshing_with_features.cpp")
create_single_source_cgal_program("test_tetrahedral_remeshing_of_one_subdomain.cpp")
create_single_source_cgal_program("test_tetrahedral_remeshing_io.cpp")
create_single_source_cgal_program("test_tetrahedral_remeshing_from_mesh_file.cpp")

# Test MLS projection
add_executable(test_tetrahedral_remeshing_mls
              "test_tetrahedral_remeshing.cpp")
target_link_libraries(test_tetrahedral_remeshing_mls PUBLIC CGAL::CGAL CGAL::Data)
target_compile_definitions(test_tetrahedral_remeshing_mls
                           PUBLIC -DCGAL_TET_REMESHING_SMOOTHING_WITH_MLS)
cgal_add_test(test_tetrahedral_remeshing_mls)


# Tests using Mesh_3 require Eigen
find_package(Eigen3 3.1.0 QUIET) #(3.1.0 or greater)
include(CGAL_Eigen3_support)
if(TARGET CGAL::Eigen3_support)
  create_single_source_cgal_program("test_mesh_and_remesh_polyhedron_with_features.cpp")
  target_link_libraries(test_mesh_and_remesh_polyhedron_with_features PUBLIC CGAL::Eigen3_support)

  create_single_source_cgal_program("test_mesh_and_remesh_with_sizing_field.cpp")
  target_link_libraries(test_mesh_and_remesh_with_sizing_field PUBLIC CGAL::Eigen3_support)

  create_single_source_cgal_program("features_and_adaptive_sizing.cpp")
  target_link_libraries(features_and_adaptive_sizing PUBLIC CGAL::Eigen3_support)

  # with MLS projection
  add_executable(test_mesh_and_remesh_polyhedron_with_features_mls
                "test_mesh_and_remesh_polyhedron_with_features.cpp")
  target_link_libraries(test_mesh_and_remesh_polyhedron_with_features_mls
                        PUBLIC CGAL::CGAL CGAL::Data CGAL::Eigen3_support)
  target_compile_definitions(test_mesh_and_remesh_polyhedron_with_features_mls
                             PUBLIC -DCGAL_TET_REMESHING_SMOOTHING_WITH_MLS)
  cgal_add_test(test_mesh_and_remesh_polyhedron_with_features_mls)

  if(CGAL_ImageIO_USE_ZLIB)
    create_single_source_cgal_program("test_mesh_and_remesh_image.cpp")
    target_link_libraries(test_mesh_and_remesh_image PUBLIC CGAL::Eigen3_support)
  else()
    message(STATUS "NOTICE: The test 'test_mesh_and_remesh_image' uses zlib, and will not be compiled.")
  endif()
else()
  message(STATUS "NOTICE: Some tests require Eigen 3.1 (or greater), and will not be compiled.")
endif()


