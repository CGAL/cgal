# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.12...3.29)
project(Vector_Graphics_on_surfaces_Examples)

# CGAL and its components
find_package(CGAL REQUIRED)

find_package(NanoSVG)
include(CGAL_NanoSVG_support)

find_package(Eigen3 3.2.0 QUIET) #(requires 3.2.0 or greater)
include(CGAL_Eigen3_support)
if(TARGET CGAL::Eigen3_support)
  create_single_source_cgal_program("trace_bezier_segment_sm_example.cpp")
  target_link_libraries(trace_bezier_segment_sm_example PUBLIC CGAL::Eigen3_support)
  create_single_source_cgal_program("geodesic_circles_sm_example.cpp")
  target_link_libraries(geodesic_circles_sm_example PUBLIC CGAL::Eigen3_support)
  create_single_source_cgal_program("straightest_geodesic_sm_example.cpp")
  target_link_libraries(straightest_geodesic_sm_example PUBLIC CGAL::Eigen3_support)
  create_single_source_cgal_program("trace_polygon_example.cpp")
  target_link_libraries(trace_polygon_example PUBLIC CGAL::Eigen3_support)
  create_single_source_cgal_program("trace_regular_polygons_example.cpp")
  target_link_libraries(trace_regular_polygons_example PUBLIC CGAL::Eigen3_support)
  create_single_source_cgal_program("trace_polygon_on_polygonal_curve_example.cpp")
  target_link_libraries(trace_polygon_on_polygonal_curve_example PUBLIC CGAL::Eigen3_support)
  if (NanoSVG_FOUND)
    create_single_source_cgal_program("trace_svg_example.cpp")
    target_link_libraries(trace_svg_example PUBLIC CGAL::Eigen3_support CGAL::NanoSVG_support)
    create_single_source_cgal_program("trace_polygon_on_svg_curve_example.cpp")
    target_link_libraries(trace_polygon_on_svg_curve_example PUBLIC CGAL::Eigen3_support CGAL::NanoSVG_support)
  else()
    message(STATUS "NOTICE: trace_svg_example requires NanoSVG and will not be compiled.")
  endif()
else()
  message(STATUS "NOTICE: Examples that use Eigen will not be compiled.")
endif()
