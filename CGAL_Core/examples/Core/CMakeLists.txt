cmake_minimum_required(VERSION 3.12...3.31)
project(Core_Examples)

# CGAL and its components
find_package(CGAL REQUIRED COMPONENTS Core)

create_single_source_cgal_program("zero-one.cpp")
create_single_source_cgal_program("delaunay.cpp")
