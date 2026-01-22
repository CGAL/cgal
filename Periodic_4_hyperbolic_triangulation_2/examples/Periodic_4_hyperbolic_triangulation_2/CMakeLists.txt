cmake_minimum_required(VERSION 3.12...3.31)
project(Periodic_4_hyperbolic_triangulation_2_Examples)

find_package(CGAL REQUIRED OPTIONAL_COMPONENTS Core)

find_package(LEDA QUIET)
create_single_source_cgal_program("p4ht2_example_insertion.cpp")

