cmake_minimum_required(VERSION 3.1...3.15)

project( Triangulation_on_sphere_2_Tests )

find_package(CGAL REQUIRED COMPONENTS Core)

find_package(Eigen3 3.1.0) #(requires 3.1.0 or greater)
include(CGAL_Eigen3_support)

if ( CGAL_FOUND )

  create_single_source_cgal_program( "test_dtos.cpp" )
  create_single_source_cgal_program( "test_dtos2_remove.cpp" )
  create_single_source_cgal_program( "test_dtos_degenerate_cases.cpp" )
  create_single_source_cgal_program( "test_dtos_illegal_points.cpp" )
  create_single_source_cgal_program( "test_dtos_projection_traits.cpp" )
  create_single_source_cgal_program( "test_dtos_traits.cpp" )

  if(TARGET CGAL::Eigen3_support)
    create_single_source_cgal_program( "test_dtos_dual.cpp" )
    target_link_libraries(test_dtos_dual PUBLIC CGAL::Eigen3_support)
  endif()

else()

    message(STATUS "This program requires the CGAL library, and will not be compiled.")

endif()
