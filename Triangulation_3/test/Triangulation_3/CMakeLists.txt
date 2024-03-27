# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.23)
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
create_single_source_cgal_program("test_RT_cell_base_with_weighted_circumcenter_3.cpp")
create_single_source_cgal_program("test_robust_weighted_circumcenter.cpp")
create_single_source_cgal_program("test_simplex_3.cpp")
create_single_source_cgal_program("test_segment_simplex_traverser_3.cpp" )
if(cxx_std_17 IN_LIST CMAKE_CXX_COMPILE_FEATURES)
  create_single_source_cgal_program("test_simplex_iterator_3.cpp" )
  create_single_source_cgal_program("test_segment_cell_traverser_3.cpp" )
  target_compile_features(test_simplex_iterator_3 PRIVATE cxx_std_17)
  target_compile_features(test_segment_simplex_traverser_3 PRIVATE cxx_std_17)
else()
  message(
  STATUS
    "NOTICE: test_simplex_iterator_3.cpp and test_segment_simplex_traverser_3.cpp require C++17 and will not be compiled."
  )
endif()
create_single_source_cgal_program("test_static_filters.cpp")
create_single_source_cgal_program("test_triangulation_3.cpp")
create_single_source_cgal_program("test_io_triangulation_3.cpp")
create_single_source_cgal_program("test_triangulation_serialization_3.cpp")
create_single_source_cgal_program("test_dt_deterministic_3.cpp")
create_single_source_cgal_program("test_Triangulation_with_transform_iterator.cpp")
create_single_source_cgal_program("test_Triangulation_with_zip_iterator.cpp")

if(TARGET CGAL::TBB_support)
  message(STATUS "Found TBB")

  foreach(target test_delaunay_3 test_regular_3
                 test_regular_insert_range_with_info)
    target_link_libraries(${target} PUBLIC CGAL::TBB_support)
  endforeach()

  if(CGAL_ENABLE_TESTING)
    set_property(TEST
      "execution   of  test_delaunay_3"
      "execution   of  test_regular_3"
      "execution   of  test_regular_insert_range_with_info"
      PROPERTY RUN_SERIAL 1)
  endif()
else()
  message(STATUS "NOTICE: The TBB library was not found. Some tests will not be available.")
endif()

if(CGAL_ENABLE_TESTING)
  set_tests_properties(
    "execution   of  test_delaunay_hierarchy_3"
    "execution   of  test_delaunay_hierarchy_3_old"
    "execution   of  test_regular_as_delaunay_3"
    "execution   of  test_triangulation_3"
    "execution   of  test_delaunay_3"
    PROPERTIES RESOURCE_LOCK Triangulation_3_Tests_IO)
endif()
