cmake_minimum_required(VERSION 3.12...3.31)
project(Constrained_triangulation_3_Tests)

find_package(CGAL REQUIRED)

create_single_source_cgal_program(test_constrained_Delaunay_triangulation_3.cpp)
create_single_source_cgal_program(test_ccdt_remeshing.cpp)
