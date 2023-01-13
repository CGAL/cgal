cmake_minimum_required(VERSION 3.1...3.23)

project(Isosurfacing_3_Tests)

find_package(CGAL REQUIRED)

create_single_source_cgal_program( "test_dual_contouring.cpp" )
create_single_source_cgal_program( "test_marching_cubes.cpp" )

find_package(TBB)
include(CGAL_TBB_support)
if(TARGET CGAL::TBB_support)
  target_link_libraries(test_dual_contouring PUBLIC CGAL::TBB_support)
  target_link_libraries(test_marching_cubes PUBLIC CGAL::TBB_support)
endif()
