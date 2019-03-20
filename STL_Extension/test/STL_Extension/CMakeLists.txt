# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( STL_Extension_Tests )

cmake_minimum_required(VERSION 3.1)


find_package(CGAL QUIET)

if ( CGAL_FOUND )

  find_package( TBB QUIET )

  create_single_source_cgal_program( "test_Boolean_tag.cpp" )
  create_single_source_cgal_program( "test_Cache.cpp" )
  create_single_source_cgal_program( "test_Compact_container.cpp" )
  create_single_source_cgal_program( "test_Compact_container_is_used.cpp" )
  create_single_source_cgal_program( "test_complexity_tags.cpp" )
  create_single_source_cgal_program( "test_composition.cpp" )
  create_single_source_cgal_program( "test_Concatenate_iterator.cpp" )
  create_single_source_cgal_program( "test_Concurrent_compact_container.cpp" )
  if(TBB_FOUND)
    CGAL_target_use_TBB(test_Concurrent_compact_container)
  endif()
  create_single_source_cgal_program( "test_dispatch_output.cpp" )
  create_single_source_cgal_program( "test_Flattening_iterator.cpp" )
  create_single_source_cgal_program( "test_Handle_with_policy.cpp" )
  create_single_source_cgal_program( "test_In_place_list.cpp" )
  create_single_source_cgal_program( "test_is_iterator.cpp" )
  create_single_source_cgal_program( "test_is_streamable.cpp" )
  create_single_source_cgal_program( "test_lexcompare_outputrange.cpp" )
  create_single_source_cgal_program( "test_Modifiable_priority_queue.cpp" )
  create_single_source_cgal_program( "test_multiset.cpp" )
  create_single_source_cgal_program( "test_N_tuple.cpp" )
  create_single_source_cgal_program( "test_Nested_iterator.cpp" )
  create_single_source_cgal_program( "test_Object.cpp" )
  create_single_source_cgal_program( "test_stl_extension.cpp" )
  create_single_source_cgal_program( "test_type_traits.cpp" )
  create_single_source_cgal_program( "test_Uncertain.cpp" )
  create_single_source_cgal_program( "test_vector.cpp" )
  create_single_source_cgal_program( "test_join_iterators.cpp" )
else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()

