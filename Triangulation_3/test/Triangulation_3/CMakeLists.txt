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
create_single_source_cgal_program("test_simplex_iterator_3.cpp" )
create_single_source_cgal_program("test_segment_cell_traverser_3.cpp" )
create_single_source_cgal_program("test_segment_simplex_traverser_3.cpp" )
if(cxx_std_17 IN_LIST CMAKE_CXX_COMPILE_FEATURES)
  target_compile_features(test_segment_simplex_traverser_3 PRIVATE cxx_std_17)
else()
  message(
  STATUS
    "NOTICE: test_segment_simplex_traverser_3.cpp requires C++17 and will not be compiled."
  )
endif()
create_single_source_cgal_program("test_static_filters.cpp")
create_single_source_cgal_program("test_triangulation_3.cpp")
create_single_source_cgal_program("test_io_triangulation_3.cpp")
create_single_source_cgal_program("test_triangulation_serialization_3.cpp")
create_single_source_cgal_program("test_dt_deterministic_3.cpp")
create_single_source_cgal_program("test_Triangulation_with_transform_iterator.cpp")
create_single_source_cgal_program("test_Triangulation_with_zip_iterator.cpp")
create_single_source_cgal_program("debug_segment_traverser_in_2D.cpp")

 # CDT_3: C++20 for structured bindings
create_single_source_cgal_program( "cdt_test_insert_constrained_edge_from_EDG_file.cpp" CXX_FEATURES cxx_std_20)
create_single_source_cgal_program( "cdt_test_insert_constrained_edge_from_OFF_file.cpp" CXX_FEATURES cxx_std_20)
create_single_source_cgal_program( "cdt_3_from_off.cpp" CXX_FEATURES cxx_std_20)
add_executable(cdt_3_from_off_CGAL_DEBUG_CDT_3 cdt_3_from_off)
target_compile_definitions(cdt_3_from_off_CGAL_DEBUG_CDT_3 PRIVATE CGAL_DEBUG_CDT_3=255)
target_compile_features(cdt_3_from_off_CGAL_DEBUG_CDT_3 PRIVATE cxx_std_20)
target_link_libraries(cdt_3_from_off_CGAL_DEBUG_CDT_3 PUBLIC CGAL::CGAL CGAL::Data)
cgal_add_test(cdt_3_from_off_CGAL_DEBUG_CDT_3)

add_executable(test_CDT_3_insert_constrained_edge_from_EDG_file cdt_test_insert_constrained_edge_from_EDG_file.cpp)
target_link_libraries(test_CDT_3_insert_constrained_edge_from_EDG_file PUBLIC CGAL::CGAL CGAL::Data)
target_compile_features(test_CDT_3_insert_constrained_edge_from_EDG_file PUBLIC cxx_std_20)
target_compile_definitions(test_CDT_3_insert_constrained_edge_from_EDG_file PUBLIC CGAL_TEST_CDT_3_USE_CDT)
cgal_add_test(test_CDT_3_insert_constrained_edge_from_EDG_file)

add_executable(test_CDT_3_insert_constrained_edge_from_OFF_file cdt_test_insert_constrained_edge_from_OFF_file.cpp)
target_link_libraries(test_CDT_3_insert_constrained_edge_from_OFF_file PUBLIC CGAL::CGAL CGAL::Data)
target_compile_features(test_CDT_3_insert_constrained_edge_from_OFF_file PUBLIC cxx_std_20)
target_compile_definitions(test_CDT_3_insert_constrained_edge_from_OFF_file PUBLIC CGAL_TEST_CDT_3_USE_CDT)
cgal_add_test(test_CDT_3_insert_constrained_edge_from_OFF_file)

add_test(NAME "execution   of  cdt_3_from_off cube"    COMMAND cdt_3_from_off ${CGAL_DATA_DIR}/meshes/cube.off    ${CMAKE_CURRENT_BINARY_DIR}/dump_output_cube.off)
add_test(NAME "execution   of  cdt_3_from_off sphere"  COMMAND cdt_3_from_off ${CGAL_DATA_DIR}/meshes/sphere.off  ${CMAKE_CURRENT_BINARY_DIR}/dump_output_sphere.off)
add_test(NAME "execution   of  cdt_3_from_off fandisk" COMMAND cdt_3_from_off ${CGAL_DATA_DIR}/meshes/fandisk.off ${CMAKE_CURRENT_BINARY_DIR}/dump_output_fandisk.off)
add_test(NAME "execution   of  cdt_3_from_off mpi"     COMMAND cdt_3_from_off ${CGAL_DATA_DIR}/meshes/mpi.off     ${CMAKE_CURRENT_BINARY_DIR}/dump_output_mpi.off)
add_test(NAME "execution   of  cdt_3_from_off 3torus"  COMMAND cdt_3_from_off ${CGAL_DATA_DIR}/meshes/3torus.off  ${CMAKE_CURRENT_BINARY_DIR}/dump_output_3torus.off)
add_test(NAME "execution   of  cdt_3_from_off_CGAL_DEBUG_CDT_3 3torus" COMMAND cdt_3_from_off_CGAL_DEBUG_CDT_3 ${CGAL_DATA_DIR}/meshes/3torus.off)
add_test(NAME "execution   of  cdt_3_from_off cheese-selection"     COMMAND cdt_3_from_off ${CGAL_DATA_DIR}/meshes/cheese-selection.off     ${CMAKE_CURRENT_BINARY_DIR}/dump_output_cheese-selection.off)
add_test(NAME "execution   of  cdt_3_from_off cheese-selection-2"     COMMAND cdt_3_from_off ${CGAL_DATA_DIR}/meshes/cheese-selection-2.off     ${CMAKE_CURRENT_BINARY_DIR}/dump_output_cheese-selection-2.off)


cgal_setup_test_properties("execution   of  cdt_3_from_off cube" cdt_3_from_off)
cgal_setup_test_properties("execution   of  cdt_3_from_off sphere" cdt_3_from_off)
cgal_setup_test_properties("execution   of  cdt_3_from_off fandisk" cdt_3_from_off)
cgal_setup_test_properties("execution   of  cdt_3_from_off mpi" cdt_3_from_off)
cgal_setup_test_properties("execution   of  cdt_3_from_off cheese-selection" cdt_3_from_off)
cgal_setup_test_properties("execution   of  cdt_3_from_off cheese-selection-2" cdt_3_from_off)
cgal_setup_test_properties("execution   of  cdt_3_from_off 3torus" cdt_3_from_off)
cgal_setup_test_properties("execution   of  cdt_3_from_off_CGAL_DEBUG_CDT_3 3torus" cdt_3_from_off_CGAL_DEBUG_CDT_3)

get_directory_property(all_tests TESTS)
foreach(test ${all_tests})
  if(test MATCHES cdt|CDT)
    set_property(TEST ${test} APPEND PROPERTY LABELS CDT_3)
  endif()
endforeach()

if(TARGET CGAL::TBB_support)
  message(STATUS "Found TBB")

  foreach(target test_delaunay_3 test_regular_3
                 test_regular_insert_range_with_info)
    target_link_libraries(${target} PUBLIC CGAL::TBB_support)
  endforeach()

  if(BUILD_TESTING)
    set_property(TEST
      "execution   of  test_delaunay_3"
      "execution   of  test_regular_3"
      "execution   of  test_regular_insert_range_with_info"
      PROPERTY RUN_SERIAL 1)
  endif()
else()
  message(STATUS "NOTICE: The TBB library was not found. Some tests will not be available.")
endif()

if(BUILD_TESTING)
  set_tests_properties(
    "execution   of  test_delaunay_hierarchy_3"
    "execution   of  test_delaunay_hierarchy_3_old"
    "execution   of  test_regular_as_delaunay_3"
    "execution   of  test_triangulation_3"
    "execution   of  test_delaunay_3"
    PROPERTIES RESOURCE_LOCK Triangulation_3_Tests_IO)
endif()
