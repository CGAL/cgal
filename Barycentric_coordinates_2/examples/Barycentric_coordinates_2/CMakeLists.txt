# Created by the script cgal_create_cmake_script.
# This is the CMake script for compiling a CGAL application.

project(Barycentric_coordinates_2_Examples)

cmake_minimum_required(VERSION 2.8.10)
set(CMAKE_CXX_STANDARD 14)

find_package(CGAL QUIET COMPONENTS Core)
if(CGAL_FOUND)

  include(${CGAL_USE_FILE})
  include(CGAL_CreateSingleSourceCGALProgram)

  create_single_source_cgal_program(
    "segment_coordinates.cpp")
  create_single_source_cgal_program(
    "triangle_coordinates.cpp")
  create_single_source_cgal_program(
    "wachspress_coordinates.cpp")
  create_single_source_cgal_program(
    "mean_value_coordinates.cpp")
  create_single_source_cgal_program(
    "discrete_harmonic_coordinates.cpp")

  # missing:
  # construct_centroid_2
  # in projection traits
  # create_single_source_cgal_program(
  #   "terrain_height_modeling.cpp")

  # this code is deprecated:
  # create_single_source_cgal_program(
  #   "deprecated_coordinates.cpp")

  find_package(Eigen3 REQUIRED)
  if(EIGEN3_FOUND)
    include(${EIGEN3_USE_FILE})

    create_single_source_cgal_program(
      "affine_coordinates.cpp")
    create_single_source_cgal_program(
      "harmonic_coordinates.cpp")
    create_single_source_cgal_program(
      "shape_deformation.cpp")

  else()
    message(WARNING
    "Harmonic/Affine coordinates require the Eigen library, and will not be compiled.")
  endif()
else()
  message(WARNING
  "This program requires the CGAL library, and will not be compiled.")
endif()
