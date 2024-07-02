# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.12...3.29)

project( Triangulation_on_sphere_2_Benchmarks )

find_package(CGAL REQUIRED COMPONENTS Core )

create_single_source_cgal_program( "bench_dtos2.cpp" )
create_single_source_cgal_program( "generate_points.cpp" )

