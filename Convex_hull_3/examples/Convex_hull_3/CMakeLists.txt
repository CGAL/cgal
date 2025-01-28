# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.12...3.29)
project(Convex_hull_3_Examples)

# CGAL and its components
find_package(CGAL REQUIRED)

create_single_source_cgal_program("quickhull_indexed_triangle_set_3.cpp")
create_single_source_cgal_program("dynamic_hull_3.cpp")
create_single_source_cgal_program("dynamic_hull_LCC_3.cpp")
create_single_source_cgal_program("dynamic_hull_SM_3.cpp")
create_single_source_cgal_program("halfspace_intersection_3.cpp")
create_single_source_cgal_program("lloyd_algorithm.cpp")
create_single_source_cgal_program("quickhull_3.cpp")
create_single_source_cgal_program("graph_hull_3.cpp")
create_single_source_cgal_program("quickhull_any_dim_3.cpp")
create_single_source_cgal_program("extreme_points_3_sm.cpp")
create_single_source_cgal_program("extreme_indices_3.cpp")

find_package(OpenMesh QUIET)
if(OpenMesh_FOUND)
  include(CGAL_OpenMesh_support)
  message(STATUS "Found OpenMesh")

  create_single_source_cgal_program("quickhull_OM_3.cpp")
  target_link_libraries(quickhull_OM_3 PRIVATE CGAL::OpenMesh_support)
  create_single_source_cgal_program("dynamic_hull_OM_3.cpp")
  target_link_libraries(dynamic_hull_OM_3 PRIVATE CGAL::OpenMesh_support)
else()
  message(STATUS "NOTICE: Examples that use OpenMesh will not be compiled.")
endif()
