cmake_minimum_required(VERSION 3.12...3.29)
project(Periodic_4_hyperbolic_triangulation_2_Examples)

find_package(CGAL REQUIRED OPTIONAL_COMPONENTS Core)

find_package(LEDA QUIET)
if(CGAL_Core_FOUND OR LEDA_FOUND)
  create_single_source_cgal_program("p4ht2_example_insertion.cpp")
else()
  message("NOTICE: This program requires the CGAL library, and will not be compiled.")
endif()
