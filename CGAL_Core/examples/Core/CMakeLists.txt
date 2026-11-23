cmake_minimum_required(VERSION 3.12...4.4)
project(Core_Examples)

# CGAL and its components
find_package(CGAL REQUIRED)

create_single_source_cgal_program("zero-one.cpp")
create_single_source_cgal_program("delaunay.cpp")
