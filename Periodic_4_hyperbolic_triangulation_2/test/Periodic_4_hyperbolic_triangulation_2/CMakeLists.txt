
cmake_minimum_required(VERSION 3.1...3.15)
project( Periodic_4_hyperbolic_triangulation_2_Tests )

find_package(CGAL REQUIRED QUIET OPTIONAL_COMPONENTS Core )
include( ${CGAL_USE_FILE} )
find_package(LEDA QUIET)

if ( CGAL_FOUND AND (CGAL_Core_FOUND OR LEDA_FOUND) )

  include( CGAL_CreateSingleSourceCGALProgram )

  foreach(cppfile
    "test_p4ht2_construct_point_2"
    "test_p4ht2_exact_complex_numbers"
    "test_p4ht2_intersections"
    "test_p4ht2_locate"
    "test_p4ht2_matrix"
    "test_p4ht2_removal_iterator"
    "test_p4ht2_removal"
    "test_p4ht2_insert_degenerate")
      create_single_source_cgal_program(${cppfile}.cpp)
      target_link_libraries(${cppfile} PRIVATE CGAL::CGAL_Core)
  endforeach()

else()

  message(STATUS "This program requires the CGAL library and the GMP library, and will not be compiled.")

endif()

