# Created by the script cgal_create_cmake_script.
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.15)
project(Shape_detection_Examples)

set(CMAKE_CXX_STANDARD 11)

find_package(CGAL QUIET COMPONENTS Core)
if(CGAL_FOUND)

  include(${CGAL_USE_FILE})
  include(CGAL_CreateSingleSourceCGALProgram)

  # Use Eigen.
  find_package(Eigen3 3.1.0 QUIET) # (3.1.0 or greater)
  if(EIGEN3_FOUND)
    include(${EIGEN3_USE_FILE})

    create_single_source_cgal_program(
      "efficient_RANSAC_basic.cpp")
    target_link_libraries(efficient_RANSAC_basic PRIVATE CGAL::CGAL_Core)
    create_single_source_cgal_program(
      "efficient_RANSAC_with_callback.cpp")
    target_link_libraries(efficient_RANSAC_with_callback PRIVATE CGAL::CGAL_Core)
    create_single_source_cgal_program(
      "efficient_RANSAC_and_plane_regularization.cpp")
    target_link_libraries(efficient_RANSAC_and_plane_regularization PRIVATE CGAL::CGAL_Core)

    create_single_source_cgal_program(
      "region_growing_on_point_set_2.cpp")
    target_link_libraries(region_growing_on_point_set_2 PRIVATE CGAL::CGAL_Core)
    create_single_source_cgal_program(
      "region_growing_on_point_set_3.cpp")
    target_link_libraries(region_growing_on_point_set_3 PRIVATE CGAL::CGAL_Core)
    create_single_source_cgal_program(
      "region_growing_on_polygon_mesh.cpp")
    target_link_libraries(region_growing_on_polygon_mesh PRIVATE CGAL::CGAL_Core)
    create_single_source_cgal_program(
      "region_growing_with_custom_classes.cpp")
    target_link_libraries(region_growing_with_custom_classes PRIVATE CGAL::CGAL_Core)

    create_single_source_cgal_program(
      "shape_detection_basic_deprecated.cpp")
    target_link_libraries(shape_detection_basic_deprecated PRIVATE CGAL::CGAL_Core)
  endif()

  create_single_source_cgal_program(
    "efficient_RANSAC_with_custom_shape.cpp")
  target_link_libraries(efficient_RANSAC_with_custom_shape PRIVATE CGAL::CGAL_Core)
  create_single_source_cgal_program(
    "efficient_RANSAC_with_parameters.cpp")
  target_link_libraries(efficient_RANSAC_with_parameters PRIVATE CGAL::CGAL_Core)
  create_single_source_cgal_program(
    "efficient_RANSAC_with_point_access.cpp")
  target_link_libraries(efficient_RANSAC_with_point_access PRIVATE CGAL::CGAL_Core)

else()
  message(WARNING
  "This program requires the CGAL library, and will not be compiled.")
endif()
