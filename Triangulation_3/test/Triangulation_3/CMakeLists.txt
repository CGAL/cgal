# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


cmake_minimum_required(VERSION 3.1...3.15)
project( Triangulation_3_Tests )



find_package(CGAL QUIET COMPONENTS Qt5)

if ( CGAL_FOUND )

  find_package( TBB QUIET )

  include_directories (BEFORE "include")

  create_single_source_cgal_program( "test_delaunay_3.cpp" )
  create_single_source_cgal_program( "test_delaunay_hierarchy_3.cpp" )
  create_single_source_cgal_program( "test_delaunay_hierarchy_3_old.cpp" )
  create_single_source_cgal_program( "test_regular_3.cpp" )
  create_single_source_cgal_program( "test_regular_as_delaunay_3.cpp" )
  create_single_source_cgal_program( "test_regular_insert_range_with_info.cpp" )
  create_single_source_cgal_program( "test_regular_remove_3.cpp" )
  create_single_source_cgal_program( "test_regular_traits_3.cpp" )
  create_single_source_cgal_program( "test_RT_cell_base_with_weighted_circumcenter_3.cpp" )
  create_single_source_cgal_program( "test_robust_weighted_circumcenter.cpp" )
  create_single_source_cgal_program( "test_simplex_3.cpp" )
  create_single_source_cgal_program( "test_simplex_iterator_3.cpp" )
  create_single_source_cgal_program( "test_segment_cell_traverser_3.cpp" )
  create_single_source_cgal_program( "test_segment_simplex_traverser_3.cpp" )
  create_single_source_cgal_program( "test_static_filters.cpp" )
  create_single_source_cgal_program( "test_triangulation_3.cpp" )

  create_single_source_cgal_program( "cdt_3.cpp" )
  target_compile_features(cdt_3 PRIVATE cxx_std_17) # Structured bindings
  if(CGAL_Qt5_FOUND)
    target_link_libraries(cdt_3 PRIVATE CGAL::CGAL_Qt5)
    target_compile_definitions(cdt_3 PRIVATE -DCGAL_USE_BASIC_VIEWER -DQT_NO_KEYWORDS)
    if(BUILD_TESTING) # @TODO Remove this special treatment of `execution___of__cdt_3`
      set_property(TEST execution___of__cdt_3 PROPERTY ENVIRONMENT)
    endif()
  endif()

  if(TBB_FOUND)
    foreach(target
        test_delaunay_3
        test_regular_3
        test_regular_insert_range_with_info)
      CGAL_target_use_TBB(${target})
    endforeach()
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

else()

    message(STATUS "This program requires the CGAL library, and will not be compiled.")

endif()

