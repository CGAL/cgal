# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


cmake_minimum_required(VERSION 3.1...3.15)
project( Periodic_3_triangulation_3_Tests )


find_package(CGAL QUIET COMPONENTS Core )

if ( CGAL_FOUND )

  include(${CGAL_USE_FILE})
  include_directories (BEFORE "include")
  foreach(cppfile
    "test_p3rt3_as_p3t3"
    "test_p3rt3_insert_remove_point_set"
    "test_p3rt3_versus_rt3"
    "test_periodic_3_alpha_shape_3"
    "test_periodic_3_delaunay_3"
    "test_periodic_3_delaunay_hierarchy_3"
    "test_periodic_3_offset_3"
    "test_periodic_3_regular_triangulation_3"
    "test_periodic_3_regular_triangulation_traits_3"
    "test_periodic_3_triangulation_3"
    "test_periodic_3_triangulation_tds"
    "test_periodic_3_triangulation_traits_C_3"
    "test_periodic_3_triangulation_traits_H_3"
    "test_periodic_3_triangulation_traits_SC_3"
    "test_periodic_3_triangulation_traits_SH_3")
      create_single_source_cgal_program(${cppfile}.cpp)
      target_link_libraries(${cppfile} PRIVATE CGAL::CGAL_Core)
  endforeach()

else()

    message(STATUS "This program requires the CGAL library, and will not be compiled.")

endif()

