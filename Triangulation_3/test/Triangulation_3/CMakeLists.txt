# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.20)
project(Triangulation_3_Tests)

find_package(CGAL REQUIRED)

find_package(TBB QUIET)
include(CGAL_TBB_support)

include_directories(BEFORE "include")

create_single_source_cgal_program("test_delaunay_3.cpp")
create_single_source_cgal_program("test_delaunay_hierarchy_3.cpp")
create_single_source_cgal_program("test_delaunay_hierarchy_3_old.cpp")
create_single_source_cgal_program("test_regular_3.cpp")
create_single_source_cgal_program("test_regular_as_delaunay_3.cpp")
create_single_source_cgal_program("test_regular_insert_range_with_info.cpp")
create_single_source_cgal_program("test_regular_remove_3.cpp")
create_single_source_cgal_program("test_regular_traits_3.cpp")
create_single_source_cgal_program(
  "test_RT_cell_base_with_weighted_circumcenter_3.cpp")
create_single_source_cgal_program("test_robust_weighted_circumcenter.cpp")
create_single_source_cgal_program("test_simplex_3.cpp")
create_single_source_cgal_program( "test_simplex_iterator_3.cpp" )
create_single_source_cgal_program( "test_segment_cell_traverser_3.cpp" )
create_single_source_cgal_program( "test_segment_simplex_traverser_3.cpp" )
create_single_source_cgal_program("test_static_filters.cpp")
create_single_source_cgal_program("test_triangulation_3.cpp")
create_single_source_cgal_program("test_io_triangulation_3.cpp")

if(TARGET CGAL::TBB_support)
  foreach(target test_delaunay_3 test_regular_3
                 test_regular_insert_range_with_info)
    target_link_libraries(${target} PUBLIC CGAL::TBB_support)
  endforeach()
  if(BUILD_TESTING)
    set_property(TEST
      execution___of__test_delaunay_3
      execution___of__test_regular_3
      execution___of__test_regular_insert_range_with_info
      PROPERTY RUN_SERIAL 1)
  endif()
endif()

if(BUILD_TESTING)
  set_tests_properties(
    execution___of__test_delaunay_hierarchy_3
    execution___of__test_delaunay_hierarchy_3_old
    execution___of__test_regular_as_delaunay_3
    execution___of__test_triangulation_3
    execution___of__test_delaunay_3
    PROPERTIES RESOURCE_LOCK Triangulation_3_Tests_IO)
endif()
