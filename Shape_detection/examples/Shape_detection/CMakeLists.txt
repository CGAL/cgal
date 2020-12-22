# Created by the script cgal_create_cmake_script.
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.15)
project(Shape_detection_Examples)

set(CMAKE_CXX_STANDARD 11)

find_package(CGAL REQUIRED COMPONENTS Core)

include(${CGAL_USE_FILE})
include(CGAL_CreateSingleSourceCGALProgram)

# Use Eigen.
find_package(Eigen3 3.1.0 QUIET) # (3.1.0 or greater)
include(CGAL_Eigen_support)
if(TARGET CGAL::Eigen_support)
  create_single_source_cgal_program("efficient_RANSAC_basic.cpp")
  create_single_source_cgal_program("efficient_RANSAC_with_callback.cpp")
  create_single_source_cgal_program(
    "efficient_RANSAC_and_plane_regularization.cpp")

  create_single_source_cgal_program("region_growing_on_point_set_2.cpp")
  create_single_source_cgal_program("region_growing_on_point_set_3.cpp")
  create_single_source_cgal_program("region_growing_on_polygon_mesh.cpp")
  create_single_source_cgal_program("region_growing_with_custom_classes.cpp")

  create_single_source_cgal_program("shape_detection_basic_deprecated.cpp")

  foreach(
    target
    efficient_RANSAC_basic
    efficient_RANSAC_with_callback
    efficient_RANSAC_and_plane_regularization
    region_growing_on_point_set_2
    region_growing_on_point_set_3
    region_growing_on_polygon_mesh
    region_growing_with_custom_classes
    shape_detection_basic_deprecated)
    target_link_libraries(${target} PUBLIC CGAL::Eigen_support)
  endforeach()

endif()

create_single_source_cgal_program("efficient_RANSAC_with_custom_shape.cpp")
create_single_source_cgal_program("efficient_RANSAC_with_parameters.cpp")
create_single_source_cgal_program("efficient_RANSAC_with_point_access.cpp")
