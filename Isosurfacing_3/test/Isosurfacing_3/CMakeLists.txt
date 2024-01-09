cmake_minimum_required(VERSION 3.1...3.23)

project(Isosurfacing_3_Tests)

find_package(CGAL REQUIRED)

# Use Eigen
find_package(Eigen3 3.2.0 REQUIRED) #(3.2.0 or greater)
include(CGAL_Eigen3_support)

# Programs
if(TARGET CGAL::Eigen3_support)
  create_single_source_cgal_program("test_dual_contouring.cpp")
  target_link_libraries(test_dual_contouring PUBLIC CGAL::Eigen3_support)
else()
  message(STATUS "test_dual_contouring.cpp requires the Eigen library, and will not be compiled.")
endif()

create_single_source_cgal_program("test_marching_cubes.cpp")

find_package(TBB)
include(CGAL_TBB_support)
if(TARGET CGAL::TBB_support)
  target_link_libraries(test_dual_contouring PUBLIC CGAL::TBB_support)
  target_link_libraries(test_marching_cubes PUBLIC CGAL::TBB_support)
endif()
