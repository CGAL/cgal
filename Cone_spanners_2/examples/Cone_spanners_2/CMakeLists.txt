cmake_minimum_required(VERSION 3.12...3.31)
project(Cone_spanners_2_Examples)

find_package(CGAL REQUIRED OPTIONAL_COMPONENTS Core)
find_package(LEDA QUIET)

create_single_source_cgal_program("compute_cones.cpp")
create_single_source_cgal_program("theta_io.cpp")
create_single_source_cgal_program("dijkstra_theta.cpp")
