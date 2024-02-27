cmake_minimum_required(VERSION 3.1...3.23)

project(Isosurfacing_3_Tests)

find_package(CGAL REQUIRED)

find_package(Eigen3 3.1.0 QUIET) #(3.1.0 or greater)
include(CGAL_Eigen3_support)

find_package(TBB QUIET)
include(CGAL_TBB_support)

create_single_source_cgal_program("test_marching_cubes.cpp")

if(TARGET CGAL::Eigen3_support)
  create_single_source_cgal_program("test_isosurfacing_concepts.cpp")
  create_single_source_cgal_program("test_dual_contouring.cpp")
  create_single_source_cgal_program("test_dual_contouring_implicit_shapes.cpp")

  target_link_libraries(test_isosurfacing_concepts PRIVATE CGAL::Eigen3_support)
  target_link_libraries(test_dual_contouring PRIVATE CGAL::Eigen3_support)
  target_link_libraries(test_dual_contouring_implicit_shapes PRIVATE CGAL::Eigen3_support)

  if(TARGET CGAL::TBB_support)
    target_link_libraries(test_isosurfacing_concepts PRIVATE CGAL::Eigen3_support)
    target_link_libraries(test_dual_contouring PRIVATE CGAL::TBB_support)
    target_link_libraries(test_dual_contouring_implicit_shapes PRIVATE CGAL::TBB_support)
  endif()
else()
  message(STATUS "NOTICE: tests require the Eigen library, and will not be compiled.")
endif()

if(TARGET CGAL::TBB_support)
  target_link_libraries(test_marching_cubes PRIVATE CGAL::TBB_support)
endif()
