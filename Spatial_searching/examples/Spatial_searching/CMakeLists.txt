# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.1...3.15)
project(Spatial_searching_Examples)

# CGAL and its components
find_package(CGAL REQUIRED)

find_package(Eigen3 3.1.91) #(requires 3.2.0 or greater)
include(CGAL_Eigen_support)

if(MSVC)
  # Turn off VC++ warning
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4244")
endif()

# include for local directory

# include for local package

# Creating entries for all .cpp/.C files with "main" routine
# ##########################################################

create_single_source_cgal_program("circular_query.cpp")

create_single_source_cgal_program("distance_browsing.cpp")

create_single_source_cgal_program("iso_rectangle_2_query.cpp")

create_single_source_cgal_program("nearest_neighbor_searching.cpp")

create_single_source_cgal_program("searching_with_circular_query.cpp")

create_single_source_cgal_program("searching_with_point_with_info.cpp")

create_single_source_cgal_program("searching_with_point_with_info_inplace.cpp")

create_single_source_cgal_program("searching_with_point_with_info_pmap.cpp")

create_single_source_cgal_program("searching_surface_mesh_vertices.cpp")

create_single_source_cgal_program("searching_polyhedron_vertices.cpp")

create_single_source_cgal_program(
  "searching_polyhedron_vertices_with_fuzzy_sphere.cpp")

create_single_source_cgal_program("user_defined_point_and_distance.cpp")

create_single_source_cgal_program("using_fair_splitting_rule.cpp")

create_single_source_cgal_program("weighted_Minkowski_distance.cpp")

if(TARGET CGAL::Eigen_support)

  create_single_source_cgal_program("fuzzy_range_query.cpp")
  target_link_libraries(fuzzy_range_query PUBLIC CGAL::Eigen_support)

  create_single_source_cgal_program("general_neighbor_searching.cpp")
  target_link_libraries(general_neighbor_searching PUBLIC CGAL::Eigen_support)

else()

  message(STATUS "fuzzy_range_query.cpp and general_neighbor_searching.cpp")
  message(
    STATUS
      "will not be compiled as they use CGAL::Epick_d which requires the Eigen library."
  )

endif()

find_package(TBB QUIET)
include(CGAL_TBB_support)
if(TARGET CGAL::TBB_support)
  create_single_source_cgal_program("parallel_kdtree.cpp")
  target_link_libraries(parallel_kdtree PUBLIC CGAL::TBB_support)
else()
  message(STATUS "parallel_kdtree.cpp requires TBB and will not be compiled")
endif()
