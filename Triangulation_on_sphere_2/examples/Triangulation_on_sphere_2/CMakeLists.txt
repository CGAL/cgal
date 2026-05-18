cmake_minimum_required(VERSION 3.12...3.31)

project( Triangulation_on_sphere_2_Examples )

find_package(CGAL REQUIRED COMPONENTS Core)
find_package(LEDA QUIET)

create_single_source_cgal_program( "triang_on_sphere.cpp" )
create_single_source_cgal_program( "triang_on_sphere_range.cpp" )
create_single_source_cgal_program( "triang_on_sphere_proj.cpp" )
create_single_source_cgal_program( "triang_on_sphere_geo.cpp" )
create_single_source_cgal_program( "triang_on_sphere_exact.cpp" )

