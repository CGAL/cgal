cmake_minimum_required(VERSION 3.1...3.23)

project(Isosurfacing_3_Tests)

find_package(CGAL REQUIRED)

find_package(Eigen3 QUIET)
include(CGAL_Eigen3_support)

find_package(TBB QUIET)
include(CGAL_TBB_support)

create_single_source_cgal_program("test_marching_cubes.cpp")
create_single_source_cgal_program("test_tmc_topology.cpp")
#create_single_source_cgal_program("verifier.cpp")
create_single_source_cgal_program("test_tmc_csg.cpp")

if(TARGET CGAL::Eigen3_support)
  create_single_source_cgal_program("test_isosurfacing_concepts.cpp")
  create_single_source_cgal_program("test_dual_contouring.cpp")

  target_link_libraries(test_isosurfacing_concepts PRIVATE CGAL::Eigen3_support)
  target_link_libraries(test_dual_contouring PRIVATE CGAL::Eigen3_support)

  if(TARGET CGAL::TBB_support)
    target_link_libraries(test_isosurfacing_concepts PRIVATE CGAL::Eigen3_support)
    target_link_libraries(test_dual_contouring PRIVATE CGAL::TBB_support)
  endif()

  #examples to be moved in example when reading to be documented
  create_single_source_cgal_program("dual_contouring_strategies.cpp")
  create_single_source_cgal_program("dual_contouring_intersection_oracles.cpp")

  target_link_libraries(dual_contouring_strategies PRIVATE CGAL::Eigen3_support)
  target_link_libraries(dual_contouring_intersection_oracles PRIVATE CGAL::Eigen3_support)

  if(TARGET CGAL::TBB_support)
    target_link_libraries(dual_contouring_strategies PRIVATE CGAL::TBB_support)
    target_link_libraries(dual_contouring_intersection_oracles PRIVATE CGAL::TBB_support)
  endif()

else()
  message(STATUS "NOTICE: tests require the Eigen library, and will not be compiled.")
endif()

if(TARGET CGAL::TBB_support)
  target_link_libraries(test_marching_cubes PRIVATE CGAL::TBB_support)
  target_link_libraries(test_tmc_topology PRIVATE CGAL::TBB_support)
  #target_link_libraries(verifier PRIVATE CGAL::TBB_support)
endif()
