# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.12...3.31)
project(PMP_Mesh_repair_Examples)

# CGAL and its components
find_package(CGAL REQUIRED)

create_single_source_cgal_program("cc_compatible_orientations.cpp")
create_single_source_cgal_program("manifoldness_repair_example.cpp")
create_single_source_cgal_program("orient_polygon_soup_example.cpp")
create_single_source_cgal_program("orientation_pipeline_example.cpp")
create_single_source_cgal_program("repair_polygon_soup_example.cpp")
create_single_source_cgal_program("stitch_borders_example.cpp" )
create_single_source_cgal_program("triangulate_polyline_example.cpp")

find_package(Eigen3 QUIET)
include(CGAL_Eigen3_support)

if(TARGET CGAL::Eigen3_support)
  create_single_source_cgal_program("hole_filling_example.cpp")
  target_link_libraries(hole_filling_example PRIVATE CGAL::Eigen3_support)
  create_single_source_cgal_program("hole_filling_example_SM.cpp")
  target_link_libraries(hole_filling_example_SM PRIVATE CGAL::Eigen3_support)
  create_single_source_cgal_program("hole_filling_visitor_example.cpp")
  target_link_libraries(hole_filling_visitor_example PRIVATE CGAL::Eigen3_support)
  create_single_source_cgal_program("hole_filling_example_LCC.cpp")
  target_link_libraries(hole_filling_example_LCC PRIVATE CGAL::Eigen3_support)
else()
  message(STATUS "NOTICE: Examples that use Eigen will not be compiled.")
endif()

find_package(OpenMesh QUIET)
if(OpenMesh_FOUND)
  include(CGAL_OpenMesh_support)

  if(TARGET CGAL::Eigen3_support)
    create_single_source_cgal_program("hole_filling_example_OM.cpp")
    target_link_libraries(hole_filling_example_OM PRIVATE CGAL::Eigen3_support CGAL::OpenMesh_support)
  endif()

  create_single_source_cgal_program("stitch_borders_example_OM.cpp")
  target_link_libraries(stitch_borders_example_OM PRIVATE CGAL::OpenMesh_support)
else()
  message(STATUS "NOTICE: Examples that use OpenMesh will not be compiled.")
endif()
