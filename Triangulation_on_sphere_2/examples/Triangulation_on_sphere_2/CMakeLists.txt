cmake_minimum_required(VERSION 3.1...3.23)

project( Triangulation_on_sphere_2_Examples )

find_package(CGAL REQUIRED COMPONENTS Core)

create_single_source_cgal_program( "triang_on_sphere.cpp" )
create_single_source_cgal_program( "triang_on_sphere_range.cpp" )
create_single_source_cgal_program( "triang_on_sphere_exact.cpp" )
create_single_source_cgal_program( "triang_on_sphere_proj.cpp" )
create_single_source_cgal_program( "triang_on_sphere_geo.cpp" )
