# Created by the script cgal_create_cmake_script.
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.22)
project(Shape_detection_Examples)

find_package(CGAL REQUIRED COMPONENTS Core)

# Use Eigen.
find_package(Eigen3 3.1.0 QUIET) # (3.1.0 or greater)
include(CGAL_Eigen3_support)
if(TARGET CGAL::Eigen3_support)
  create_single_source_cgal_program("efficient_RANSAC_basic.cpp")
  create_single_source_cgal_program("efficient_RANSAC_with_callback.cpp")

  create_single_source_cgal_program("region_growing_lines_on_point_set_2.cpp")
  create_single_source_cgal_program("region_growing_planes_on_point_set_3.cpp")
  create_single_source_cgal_program("region_growing_cylinders_on_point_set_3.cpp")
  create_single_source_cgal_program("region_growing_spheres_on_point_set_3.cpp")
  create_single_source_cgal_program("region_growing_circles_on_point_set_2.cpp")
  create_single_source_cgal_program("region_growing_planes_on_polygon_mesh.cpp")
  create_single_source_cgal_program("region_growing_with_custom_classes.cpp")

  create_single_source_cgal_program("shape_detection_basic_deprecated.cpp")

  foreach(
    target
    efficient_RANSAC_basic
    efficient_RANSAC_with_callback
    region_growing_lines_on_point_set_2
    region_growing_planes_on_point_set_3
    region_growing_cylinders_on_point_set_3
    region_growing_spheres_on_point_set_3
    region_growing_circles_on_point_set_2
    region_growing_planes_on_polygon_mesh
    region_growing_with_custom_classes
    shape_detection_basic_deprecated)
    target_link_libraries(${target} PUBLIC CGAL::Eigen3_support)
  endforeach()

endif()

create_single_source_cgal_program("efficient_RANSAC_with_custom_shape.cpp")
create_single_source_cgal_program("efficient_RANSAC_with_parameters.cpp")
create_single_source_cgal_program("efficient_RANSAC_with_point_access.cpp")
