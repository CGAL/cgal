cmake_minimum_required(VERSION 3.1...3.15)
project(Periodic_4_hyperbolic_triangulation_2_Examples)

find_package(CGAL REQUIRED QUIET OPTIONAL_COMPONENTS Core)
include(${CGAL_USE_FILE})
find_package(LEDA QUIET)

if((CGAL_Core_FOUND OR LEDA_FOUND))

  include(CGAL_CreateSingleSourceCGALProgram)

  create_single_source_cgal_program("p4ht2_example_insertion.cpp")

else()

  message(
    STATUS "This program requires the CGAL library, and will not be compiled.")

endif()
