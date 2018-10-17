# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Triangulation_2_Tests )

cmake_minimum_required(VERSION 2.8.10)

find_package(CGAL QUIET)

if ( CGAL_FOUND )

  include_directories (BEFORE "include")

  # create a target per cppfile
  file(GLOB cppfiles RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
  foreach(cppfile ${cppfiles})
    create_single_source_cgal_program( "${cppfile}" )
  endforeach()

  if(BUILD_TESTING)
    set_tests_properties(
        execution___of__test_constrained_triangulation_2
        execution___of__test_delaunay_triangulation_2
        execution___of__test_triangulation_geom_traits
        execution___of__test_triangulation_2
        execution___of__test_triangulation_2_bis
        execution___of__test_delaunay_hierarchy_2
        execution___of__test_const_triang_plus_2
        execution___of__test_regular_triangulation_2
        execution___of__test_const_del_triangulation_2
        execution___of__test_regular_hierarchy_2
      PROPERTIES RESOURCE_LOCK Triangulation_2_Tests_IO)
  endif()

else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()

