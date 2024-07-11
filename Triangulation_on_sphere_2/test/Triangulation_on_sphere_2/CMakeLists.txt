cmake_minimum_required(VERSION 3.12...3.29)

project( Triangulation_on_sphere_2_Tests )

find_package(CGAL REQUIRED COMPONENTS Core)
find_package(LEDA QUIET)

create_single_source_cgal_program( "test_dtos.cpp" )
create_single_source_cgal_program( "test_dtos2_remove.cpp" )
create_single_source_cgal_program( "test_dtos_illegal_points.cpp" )
create_single_source_cgal_program( "test_dtos_projection_traits.cpp" )
create_single_source_cgal_program( "test_dtos_traits.cpp" )
create_single_source_cgal_program( "issue_8200.cpp" )

if (CGAL_Core_FOUND OR LEDA_FOUND)
  create_single_source_cgal_program( "test_dtos_degenerate_cases.cpp" )
else()
  message("NOTICE: Some tests require CGAL_Core (or LEDA), and will not be compiled.")
endif()

find_package(Eigen3 3.1.0 QUIET) #(requires 3.1.0 or greater)
include(CGAL_Eigen3_support)
if(TARGET CGAL::Eigen3_support)
  create_single_source_cgal_program( "test_dtos_dual.cpp" )
  target_link_libraries(test_dtos_dual PUBLIC CGAL::Eigen3_support)
else()
  message(STATUS "NOTICE: The Eigen library was not found. The test 'test_dtos_dual' will not be compiled.")
endif()
