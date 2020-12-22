cmake_minimum_required(VERSION 3.1...3.15)
project(Optimal_transportation_reconstruction_2_Examples)

find_package(CGAL REQUIRED)

create_single_source_cgal_program("otr2_simplest_example.cpp")
create_single_source_cgal_program("otr2_simplest_example_with_tolerance.cpp")
create_single_source_cgal_program("otr2_list_output_example.cpp")
create_single_source_cgal_program("otr2_indexed_output_example.cpp")
create_single_source_cgal_program("otr2_mass_example.cpp")
