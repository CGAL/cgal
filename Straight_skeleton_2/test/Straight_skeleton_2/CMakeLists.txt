# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.23)
project(Straight_skeleton_2_Tests)

find_package(CGAL REQUIRED COMPONENTS Qt5 Core)
find_package(LEDA QUIET)

include_directories(BEFORE "include")

create_single_source_cgal_program("issue4533.cpp")
create_single_source_cgal_program("issue4684.cpp")
create_single_source_cgal_program("test_sls.cpp")
create_single_source_cgal_program("test_sls_previous_issues.cpp")
create_single_source_cgal_program("test_sls_traits.cpp")
create_single_source_cgal_program("test_straight_skeleton_copy.cpp")

if(CGAL_Qt5_FOUND)
  target_link_libraries(issue4684 PRIVATE CGAL::CGAL_Basic_viewer)
  target_link_libraries(test_sls_previous_issues PRIVATE CGAL::CGAL_Basic_viewer)
endif()

if (CGAL_Core_FOUND OR LEDA_FOUND)
  create_single_source_cgal_program("issue7149.cpp")
  create_single_source_cgal_program("issue7284.cpp")
  create_single_source_cgal_program("test_sls_offset.cpp")
  create_single_source_cgal_program("test_sls_weighted_offset.cpp")
  create_single_source_cgal_program("test_sls_simple.cpp")
  create_single_source_cgal_program("test_sls_weighted_polygons.cpp")
  create_single_source_cgal_program("test_sls_weighted_polygons_with_holes.cpp")
  if(CGAL_Qt5_FOUND)
    target_link_libraries(issue7149 PRIVATE CGAL::CGAL_Basic_viewer)
    target_link_libraries(issue7284 PRIVATE CGAL::CGAL_Basic_viewer)
    target_link_libraries(test_sls_offset PRIVATE CGAL::CGAL_Basic_viewer)
    target_link_libraries(test_sls_weighted_offset PRIVATE CGAL::CGAL_Basic_viewer)
    target_link_libraries(test_sls_weighted_polygons PRIVATE CGAL::CGAL_Basic_viewer)
    target_link_libraries(test_sls_weighted_polygons_with_holes PRIVATE CGAL::CGAL_Basic_viewer)
  endif()
endif()