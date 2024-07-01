cmake_minimum_required(VERSION 3.1...3.23)
project(Constrained_triangulation_3_Examples)

find_package(CGAL REQUIRED COMPONENTS Qt6)

create_single_source_cgal_program(constrained_Delaunay_triangulation_3.cpp)
create_single_source_cgal_program(remesh_constrained_Delaunay_triangulation_3.cpp)

target_compile_features(constrained_Delaunay_triangulation_3 PUBLIC cxx_std_20) # TODO: remove this line
target_compile_features(remesh_constrained_Delaunay_triangulation_3 PUBLIC cxx_std_20)

if(CGAL_Qt6_FOUND)
  target_link_libraries(constrained_Delaunay_triangulation_3 PUBLIC CGAL::CGAL_Basic_viewer)
else()
  message(STATUS "NOTICE: The example 'constrained_delaunay_triangulation_3' requires Qt6, and will not be compiled.")
endif()
