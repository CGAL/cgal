# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.12...3.31)
project(Straight_skeleton_extrusion_2_Tests)

find_package(CGAL REQUIRED COMPONENTS Qt6 Core)
find_package(LEDA QUIET)

include_directories(BEFORE "include")

create_single_source_cgal_program("test_sls_extrude.cpp")
create_single_source_cgal_program("test_sls_extrude_kernel_choice.cpp")
create_single_source_cgal_program("test_sls_extrude_kernels.cpp")
create_single_source_cgal_program("test_sls_extrude_degenerate.cpp")
create_single_source_cgal_program("test_sls_extrude_inward.cpp")
create_single_source_cgal_program("test_sls_extrude_outward.cpp")
create_single_source_cgal_program("test_sls_extrude_maximum_height.cpp")
create_single_source_cgal_program("test_sls_extrude_invalid_inputs.cpp")
create_single_source_cgal_program("test_sls_extrude_random.cpp")
create_single_source_cgal_program("test_sls_extrude_weighted_polygons_with_holes.cpp")

create_single_source_cgal_program("test_sls_extrude_snapping.cpp")

if(CGAL_Qt6_FOUND)
  target_link_libraries(test_sls_extrude PRIVATE CGAL::CGAL_Basic_viewer)
  target_link_libraries(test_sls_extrude_degenerate PRIVATE CGAL::CGAL_Basic_viewer)
  target_link_libraries(test_sls_extrude_invalid_inputs PRIVATE CGAL::CGAL_Basic_viewer)
  target_link_libraries(test_sls_extrude_inward PRIVATE CGAL::CGAL_Basic_viewer)
  target_link_libraries(test_sls_extrude_outward PRIVATE CGAL::CGAL_Basic_viewer)
  target_link_libraries(test_sls_extrude_kernels PRIVATE CGAL::CGAL_Basic_viewer)
  target_link_libraries(test_sls_extrude_maximum_height PRIVATE CGAL::CGAL_Basic_viewer)
  target_link_libraries(test_sls_extrude_random PRIVATE CGAL::CGAL_Basic_viewer)
  target_link_libraries(test_sls_extrude_snapping PRIVATE CGAL::CGAL_Basic_viewer)
  target_link_libraries(test_sls_extrude_weighted_polygons_with_holes PRIVATE CGAL::CGAL_Basic_viewer)
endif()
