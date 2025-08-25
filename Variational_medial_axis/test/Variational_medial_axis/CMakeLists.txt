cmake_minimum_required(VERSION 3.12...3.31)

project( Variational_medial_axis_Tests )

find_package(CGAL REQUIRED)
find_package(Eigen3 REQUIRED 3.2.0)

include(CGAL_Eigen3_support)
if(TARGET CGAL::Eigen3_support)
  create_single_source_cgal_program("variational_medial_axis_sampling_test.cpp")
  target_link_libraries(variational_medial_axis_sampling_test PRIVATE CGAL::Eigen3_support)
else()
  message("NOTICE: These tests require the Eigen library (3.2 or greater), and will not be compiled.")
endif()
