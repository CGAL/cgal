cmake_minimum_required(VERSION 3.1...3.23)
project(Triangulation_3_Examples)

find_package(CGAL REQUIRED OPTIONAL_COMPONENTS Qt6)

include_directories(BEFORE "../../include")

create_single_source_cgal_program("for_loop.cpp")
create_single_source_cgal_program("adding_handles_3.cpp")
create_single_source_cgal_program("color.cpp")
create_single_source_cgal_program("copy_triangulation_3.cpp")
create_single_source_cgal_program("fast_location_3.cpp")
create_single_source_cgal_program("find_conflicts_3.cpp")
create_single_source_cgal_program("info_insert_with_pair_iterator.cpp")
create_single_source_cgal_program("info_insert_with_pair_iterator_regular.cpp")
create_single_source_cgal_program("info_insert_with_transform_iterator.cpp")
create_single_source_cgal_program("info_insert_with_zip_iterator.cpp")
create_single_source_cgal_program("regular_3.cpp")
create_single_source_cgal_program( "segment_cell_traverser_3.cpp" )
create_single_source_cgal_program( "segment_simplex_traverser_3.cpp" )
create_single_source_cgal_program("regular_with_info_3.cpp")
create_single_source_cgal_program("simple_triangulation_3.cpp")
create_single_source_cgal_program("simplex.cpp")

create_single_source_cgal_program("draw_triangulation_3.cpp")
if(CGAL_Qt6_FOUND)
  target_link_libraries(draw_triangulation_3 PUBLIC CGAL::CGAL_Basic_viewer)
else()
  message(STATUS "NOTICE: The example 'draw_triangulation_3' requires Qt6, and will not be compiled.")
endif()

find_package(TBB QUIET)
include(CGAL_TBB_support)

if(TARGET CGAL::TBB_support)
  message(STATUS "Found TBB")

  create_single_source_cgal_program("parallel_insertion_and_removal_in_regular_3.cpp")
  create_single_source_cgal_program("parallel_insertion_in_delaunay_3.cpp")
  create_single_source_cgal_program("sequential_parallel.cpp")
  target_link_libraries(parallel_insertion_and_removal_in_regular_3 PUBLIC CGAL::TBB_support)
  target_link_libraries(parallel_insertion_in_delaunay_3 PUBLIC CGAL::TBB_support)
  target_link_libraries(sequential_parallel PUBLIC CGAL::TBB_support)

  if(CGAL_ENABLE_TESTING)
    set_property(TEST
      "execution   of  parallel_insertion_and_removal_in_regular_3"
      "execution   of  parallel_insertion_in_delaunay_3"
      "execution   of  sequential_parallel"
      PROPERTY RUN_SERIAL 1)
  endif()
else()
  message(STATUS "NOTICE: A few examples require TBB, and will not be compiled.")
endif()
