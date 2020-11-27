# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.1...3.15)
project(Point_set_3_Examples)

# CGAL and its components
find_package(CGAL REQUIRED OPTIONAL_COMPONENTS Qt5)

if(CGAL_Qt5_FOUND)
  add_definitions(-DCGAL_USE_BASIC_VIEWER -DQT_NO_KEYWORDS)
endif()

# Boost and its components
find_package(Boost REQUIRED)

if(NOT Boost_FOUND)

  message(
    STATUS "This project requires the Boost library, and will not be compiled.")

  return()

endif()

# Creating entries for all C++ files with "main" routine
# ##########################################################

create_single_source_cgal_program("point_set.cpp")
create_single_source_cgal_program("point_set_property.cpp")
create_single_source_cgal_program("point_set_read_xyz.cpp")
create_single_source_cgal_program("point_set_advanced.cpp")

set(needed_cxx_features cxx_rvalue_references cxx_variadic_templates)
create_single_source_cgal_program("point_set_read_ply.cpp" CXX_FEATURES
                                  ${needed_cxx_features})

find_package(Eigen3 3.1.0) #(requires 3.1.0 or greater)
include(CGAL_Eigen_support)
if(EIGEN3_FOUND)
  create_single_source_cgal_program("point_set_algo.cpp")
  target_link_libraries(point_set_algo PUBLIC CGAL::Eigen_support)
endif()

create_single_source_cgal_program("draw_point_set_3.cpp")
if(CGAL_Qt5_FOUND)
  target_link_libraries(draw_point_set_3 PUBLIC CGAL::CGAL_Qt5)
endif()
