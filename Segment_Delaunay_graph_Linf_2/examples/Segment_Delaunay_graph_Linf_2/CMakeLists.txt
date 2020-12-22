# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.15)
project(Segment_Delaunay_graph_Linf_2_Examples)

find_package(CGAL REQUIRED)

create_single_source_cgal_program("sdg-count-sites-linf.cpp")
create_single_source_cgal_program("sdg-fast-sp-linf.cpp")
create_single_source_cgal_program("sdg-fast-sp-polygon-linf.cpp")
create_single_source_cgal_program("sdg-filtered-traits-linf.cpp")
create_single_source_cgal_program("sdg-info-set-linf.cpp")
create_single_source_cgal_program("sdg-red-blue-info-linf.cpp")
create_single_source_cgal_program("sdg-voronoi-edges-linf.cpp")
create_single_source_cgal_program("sdg-voronoi-edges-exact-linf.cpp")
create_single_source_cgal_program("print-sdg-linf.cpp")
