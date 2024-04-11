cmake_minimum_required(VERSION 3.1...3.23)
project(Cone_spanners_2_Examples)

find_package(CGAL REQUIRED OPTIONAL_COMPONENTS Core)
find_package(LEDA QUIET)

if(CGAL_Core_FOUND OR LEDA_FOUND)
  create_single_source_cgal_program("compute_cones.cpp")
  create_single_source_cgal_program("theta_io.cpp")
else()
  message("NOTICE: This program requires the CGAL_Core library (or LEDA), and will not be compiled.")
endif()

create_single_source_cgal_program("dijkstra_theta.cpp")
