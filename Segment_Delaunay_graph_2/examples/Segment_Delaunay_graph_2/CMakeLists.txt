# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.23)
project(Segment_Delaunay_graph_2_Examples)

find_package(CGAL REQUIRED OPTIONAL_COMPONENTS Core)

create_single_source_cgal_program("sdg-count-sites.cpp")
create_single_source_cgal_program("sdg-fast-sp.cpp")
create_single_source_cgal_program("sdg-fast-sp-polygon.cpp")
create_single_source_cgal_program("sdg-filtered-traits.cpp")
create_single_source_cgal_program("sdg-info-set.cpp")
create_single_source_cgal_program("sdg-red-blue-info.cpp")
create_single_source_cgal_program("sdg-voronoi-edges.cpp")

find_package(LEDA QUIET)
if(CGAL_Core_FOUND OR LEDA_FOUND)
  create_single_source_cgal_program("sdg-advanced-draw.cpp")
else()
  message("NOTICE: The program sdg-advanced-draw requires CGAL_Core or LEDA, and will not be compiled.")
endif()
