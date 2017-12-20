project( Triangulation_3_Examples )

cmake_minimum_required(VERSION 2.8.12)


find_package(CGAL QUIET)

if ( CGAL_FOUND )
  create_single_source_cgal_program( "adding_handles_3.cpp" )
  create_single_source_cgal_program( "color.cpp" )
  create_single_source_cgal_program( "copy_triangulation_3.cpp" )
  create_single_source_cgal_program( "fast_location_3.cpp" )
  create_single_source_cgal_program( "find_conflicts_3.cpp" )
  create_single_source_cgal_program( "info_insert_with_pair_iterator.cpp" )
  create_single_source_cgal_program( "info_insert_with_pair_iterator_regular.cpp" )
  create_single_source_cgal_program( "info_insert_with_transform_iterator.cpp" )
  create_single_source_cgal_program( "info_insert_with_zip_iterator.cpp" )
  create_single_source_cgal_program( "regular_3.cpp" )
  create_single_source_cgal_program( "regular_with_info_3.cpp" )
  create_single_source_cgal_program( "simple_triangulation_3.cpp" )
  create_single_source_cgal_program( "simplex.cpp" )

  find_package( TBB QUIET )

  if( TBB_FOUND )
    include( CGAL_target_use_TBB )

    create_single_source_cgal_program( "parallel_insertion_and_removal_in_regular_3.cpp" )
    create_single_source_cgal_program( "parallel_insertion_in_delaunay_3.cpp" )
    create_single_source_cgal_program( "sequential_parallel.cpp" )
    CGAL_target_use_TBB( parallel_insertion_and_removal_in_regular_3 )
    CGAL_target_use_TBB( parallel_insertion_in_delaunay_3 )
    CGAL_target_use_TBB( sequential_parallel )
  else()
    message(STATUS "NOTICE: a few examples require TBB and will not be compiled.")
  endif()

else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()

