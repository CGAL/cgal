cmake_minimum_required(VERSION 3.12...3.31)
project(Mesh_smoothing_3_Tests)

find_package(CGAL REQUIRED)
find_package(TBB QUIET)
find_package(Eigen3 QUIET)
include(CGAL_Eigen3_support)
find_package(OpenMP COMPONENTS CXX)

if(NOT TARGET CGAL::Eigen3_support)
  message(FATAL_ERROR "Mesh_smoothing_3 tests require Eigen.")
endif()

create_single_source_cgal_program("test_boundary_aware_mesh_smoothing.cpp")
create_single_source_cgal_program("test_boundary_aware_mesh_smoothing_projectors.cpp")
create_single_source_cgal_program("test_boundary_aware_mesh_smoothing_integration.cpp")

target_link_libraries(test_boundary_aware_mesh_smoothing PRIVATE CGAL::Eigen3_support)
target_link_libraries(test_boundary_aware_mesh_smoothing_projectors PRIVATE CGAL::Eigen3_support)
target_link_libraries(test_boundary_aware_mesh_smoothing_integration PRIVATE CGAL::Eigen3_support)
target_compile_definitions(test_boundary_aware_mesh_smoothing_integration PRIVATE
  MS3_EXAMPLE_MESH_FILE="${CMAKE_CURRENT_SOURCE_DIR}/data/sphere.mesh")


if(TARGET CGAL::TBB_support)
  target_link_libraries(test_boundary_aware_mesh_smoothing PRIVATE CGAL::TBB_support)
  target_link_libraries(test_boundary_aware_mesh_smoothing_projectors PRIVATE CGAL::TBB_support)
  target_link_libraries(test_boundary_aware_mesh_smoothing_integration PRIVATE CGAL::TBB_support)
endif()

if(OpenMP_CXX_FOUND)
  if(MSVC)
      # Override the compile options for the imported OpenMP target
      set_target_properties(OpenMP::OpenMP_CXX PROPERTIES
                          INTERFACE_COMPILE_OPTIONS "/openmp:llvm")
  endif()
  target_link_libraries(test_boundary_aware_mesh_smoothing_integration PRIVATE OpenMP::OpenMP_CXX)
  target_link_libraries(test_boundary_aware_mesh_smoothing_projectors PRIVATE OpenMP::OpenMP_CXX)
  target_link_libraries(test_boundary_aware_mesh_smoothing PRIVATE OpenMP::OpenMP_CXX)
endif()
