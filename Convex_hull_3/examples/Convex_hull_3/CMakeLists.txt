# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.1...3.15)
project(Convex_hull_3_Examples)

# CGAL and its components
find_package(CGAL REQUIRED)

# Boost and its components
find_package(Boost REQUIRED)

if(NOT Boost_FOUND)

  message(
    STATUS "This project requires the Boost library, and will not be compiled.")

  return()

endif()

find_package(OpenMesh QUIET)

if(OpenMesh_FOUND)
  include(UseOpenMesh)
else()
  message(STATUS "Examples that use OpenMesh will not be compiled.")
endif()

# include for local directory

# include for local package

# Creating entries for all C++ files with "main" routine
# ##########################################################

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

if(OpenMesh_FOUND)
  create_single_source_cgal_program("quickhull_OM_3.cpp")
  create_single_source_cgal_program("dynamic_hull_OM_3.cpp")

  target_link_libraries(quickhull_OM_3 PRIVATE ${OPENMESH_LIBRARIES})
  target_link_libraries(dynamic_hull_OM_3 PRIVATE ${OPENMESH_LIBRARIES})
endif()
